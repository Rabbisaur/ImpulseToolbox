# NSx Fast Readers

This directory contains optimized NSx file readers for both MATLAB and Python.

## Performance Comparison

For a 9.9M sample file with 144 channels:

| Implementation | Load Time | Speedup vs openNSx |
|----------------|-----------|-------------------|
| **MATLAB openNSx** (baseline) | 94.5s | 1x |
| **Python open_nsx_fast** | 8.1s | **~12x** |
| **MATLAB openNSxFast** | 0.58s | **~163x** |

## MATLAB: openNSxFast.m

**Location**: `C:\Users\wang\Documents\GitHub\ImpulseToolbox\cerebus\NPMK\openNSxFast.m`

### Usage

```matlab
% Load with UV conversion
NSx = openNSxFast('filename.ns6', 'uv');

% Load raw data
NSx = openNSxFast('filename.ns6');
```

### Features

- **163x faster than openNSx**
- Exact output match with openNSx (validated)
- Supports NEURALCD (specs 2.2/2.3/3.0) and NEURALSG (spec 2.1)
- Single vectorized `fread` operation
- Minimal memory overhead

## Python: open_nsx_fast.py

**Location**: `C:\Users\wang\Documents\GitHub\ImpulseToolbox\cerebus\NPMK\open_nsx_fast.py`

### Usage

```python
from open_nsx_fast import open_nsx_fast

# Load with UV conversion
nsx = open_nsx_fast('filename.ns6', convert_to_uv=True)

# Access data
data = nsx.data  # NumPy array [channels x samples]
fs = nsx.meta_tags.sampling_freq
channel_ids = nsx.meta_tags.channel_ids
```

### Features

- **~12x faster than MATLAB openNSx**
- Pure Python with NumPy (no compilation needed)
- Returns structured data with dataclasses
- Same optimization principles as MATLAB version
- Supports both file formats

## Key Optimizations

Both implementations achieve performance through:

1. **Vectorized I/O**: Single bulk read operation instead of chunked reads
2. **Pre-allocated arrays**: No dynamic resizing during load
3. **Minimal data copying**: Direct type conversion without intermediate buffers
4. **Efficient header parsing**: Batch read + struct unpacking instead of multiple seeks

## When to Use Which

- **MATLAB openNSxFast**: Maximum performance for MATLAB workflows
- **Python open_nsx_fast**: Python data pipelines, good performance without MEX
- **MATLAB openNSx**: Compatibility with complex features (time ranges, channel selection)
