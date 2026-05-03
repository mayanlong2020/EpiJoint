import torch
import functools
_orig_load = torch.load
@functools.wraps(_orig_load)
def _patched_load(*args, **kwargs):
    if 'weights_only' in kwargs: kwargs['weights_only'] = False
    return _orig_load(*args, **kwargs)
torch.load = _patched_load
print('  [Conda-Patch] Torch.load bypass active.')
