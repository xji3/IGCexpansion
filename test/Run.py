from IGCexpansion.CodonGeneconv import ReCodonGeneconv

if __name__ == '__main__':
    paralog = ['YLR406C', 'YDL075W']
    Force = None
    alignment_file = './YLR406C_YDL075W_input.fasta'
    newicktree = './YeastTree.newick'
#     if args.force:
#         if args.model == 'MG94':
#             Force = {5:0.0}
#         elif args.model == 'HKY':
#             Force = {4:0.0}
#     else:
#         Force = None
    Force = None
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = None)
    test.get_mle(True, True, 0, 'BFGS')
    #test.get_individual_summary(summary_path = './Summary/')
    #test.get_SitewisePosteriorSummary(summary_path = './Summary/')
