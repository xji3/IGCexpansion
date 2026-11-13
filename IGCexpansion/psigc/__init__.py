"""psigc: fast composite likelihood and exact gradients for the HKY + PS-IGC pair-site model."""
from .tree import Tree
from .data import Dataset
from .model import Params, unpack_x, pack_x, hky_Q
from .engine import PSEngine, ISEngine

__all__ = ['Tree', 'Dataset', 'Params', 'unpack_x', 'pack_x', 'hky_Q', 'PSEngine', 'ISEngine']
__version__ = '0.2.0'
