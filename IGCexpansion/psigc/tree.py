"""Rooted tree with one duplication node.

Nodes in the subtree of the duplication node (including it) carry two paralogs ("post-duplication");
the other nodes carry one copy ("pre-duplication"). Edge order and the legacy x_rate layout follow
IGCexpansion's Tree / PSJSGeneconv: edges are sorted (parent, child) name tuples; if the root is not
the duplication node, the two root edges share one parameter, log of their total length, split
0.1 : 0.9 between root -> duplication node and root -> outgroup.
"""
import numpy as np


def parse_newick(s):
    """Minimal Newick parser: returns (root_name, children dict). Branch lengths are ignored."""
    s = s.strip()
    if ';' in s:
        s = s[:s.index(';')]
    else:                                    # legacy files: "(...)N0   0"
        k = s.rindex(')')
        s = s[:k + 1] + s[k + 1:].split()[0] if s[k + 1:].split() else s[:k + 1]
    s = ''.join(s.split())
    pos = 0
    children = {}
    counter = [0]

    def name_at():
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in ',();':
            pos += 1
        tok = s[start:pos].strip()
        return tok.split(':')[0]

    def node():
        nonlocal pos
        kids = []
        if s[pos] == '(':
            pos += 1
            while True:
                kids.append(node())
                if s[pos] == ',':
                    pos += 1
                    continue
                if s[pos] == ')':
                    pos += 1
                    break
        nm = name_at()
        if not nm:
            counter[0] += 1
            nm = f'_anon{counter[0]}'
        children[nm] = kids
        return nm

    root = node()
    return root, children


class Tree:
    def __init__(self, newick, dup_edge=None, dup_name='D1', tie_root=True):
        """dup_edge: (parent, child) on which the duplication node is inserted (default: root ->
        the root's first internal child). With dup_edge = None and the root itself named dup_name,
        the root is the duplication node."""
        if newick.strip().startswith('(') is False:
            newick = open(newick).read()
        root, ch = parse_newick(newick)
        if root != dup_name:
            if dup_edge is None:
                internal = [c for c in ch[root] if ch[c]]
                dup_edge = (root, internal[0])
            a, b = dup_edge
            ch[a] = [dup_name if c == b else c for c in ch[a]]
            ch[dup_name] = [b]
        self.root, self.children, self.dup = root, ch, dup_name
        self.parent = {c: p for p, cs in ch.items() for c in cs}
        self.tips = [v for v in ch if not ch[v]]
        post = set()
        stack = [dup_name]
        while stack:
            v = stack.pop(); post.add(v); stack.extend(ch[v])
        self.is_post = {v: v in post for v in ch}
        self.post_tips = [v for v in self.tips if self.is_post[v]]
        self.pre_tips = [v for v in self.tips if not self.is_post[v]]
        order = []

        def visit(v):
            for c in ch[v]:
                visit(c)
            order.append(v)
        visit(root)
        self.postorder = order
        self.edges = sorted((p, c) for p, cs in ch.items() for c in cs)
        self.post_edges = [e for e in self.edges if self.is_post[e[0]]]
        self.pre_edges = [e for e in self.edges if not self.is_post[e[0]]]
        self.root_by_dup = root == dup_name
        self.tie_root = tie_root and not self.root_by_dup
        if self.tie_root:
            re = [e for e in self.edges if e[0] == root]
            assert len(re) == 2 and any(e[1] == dup_name for e in re), 'tied root needs exactly root->D1 and one outgroup edge'
            self.free_edges = [e for e in self.edges if e[0] != root]
        else:
            self.free_edges = list(self.edges)
        self.n_rate = len(self.free_edges) + (1 if self.tie_root else 0)

    # ---- legacy x_rate layout
    def blen_from_x(self, xr):
        xr = np.asarray(xr, float)
        assert len(xr) == self.n_rate
        blen = {}
        if self.tie_root:
            tot = np.exp(xr[0])
            for e in self.edges:
                if e[0] == self.root:
                    blen[e] = tot * (0.1 if e[1] == self.dup else 0.9)
            rest = xr[1:]
        else:
            rest = xr
        for e, v in zip(self.free_edges, rest):
            blen[e] = float(np.exp(v))
        return blen

    def x_from_blen(self, blen):
        out = []
        if self.tie_root:
            out.append(np.log(sum(blen[e] for e in self.edges if e[0] == self.root)))
        out.extend(np.log(blen[e]) for e in self.free_edges)
        return np.array(out)

    def dlogblen_to_x(self, g_edge):
        """Map d lnL / d log t_e (dict) to the x_rate gradient (tied root edges add)."""
        out = []
        if self.tie_root:
            out.append(sum(g_edge.get(e, 0.0) for e in self.edges if e[0] == self.root))
        out.extend(g_edge.get(e, 0.0) for e in self.free_edges)
        return np.array(out)
