# Audit Report: LoadMatFast.m & SaveMatFast.m

**Date**: 2026-02-12
**Files Audited**: LoadMatFast.m, SaveMatFast.m
**Purpose**: Fast binary file I/O for MATLAB data

---

## Executive Summary

These files implement a custom binary format for faster save/load operations compared to MATLAB's built-in `save()` function. The implementation is functional but has significant limitations that prevent it from being a drop-in replacement for MATLAB's workspace save/load.

**Grade**: C+ (Functional but limited)

---

## Current Capabilities

### ✅ **Strengths**

1. **Fast Binary I/O** - Direct binary read/write without MATLAB's overhead
2. **Recursive Serialization** - Supports complex nested structures
3. **Type Support** - Handles 12 numeric types + struct + cell
4. **Header Metadata** - Stores dimensions, type, date, comment
5. **Error Checking** - Checksum validation for simple types

### Supported Data Types

| Type | Supported | Notes |
|------|-----------|-------|
| double, single | ✅ | Full support |
| int8-64, uint8-64 | ✅ | All integer types |
| char, logical | ✅ | Basic types |
| struct | ✅ | Recursive serialization |
| cell | ✅ | Recursive serialization |
| **Multiple variables** | ❌ | **Major limitation** |

---

## Critical Limitations

### ❌ **1. Single Variable Only**

**Current**:
```matlab
SaveMatFast('data.bin', myVariable, 'comment')
data = LoadMatFast('data.bin');
```

**MATLAB standard**:
```matlab
save('data.mat', 'var1', 'var2', 'var3')
load('data.mat')  % loads all variables into workspace
```

**Impact**: Cannot replace MATLAB's save/load for multi-variable workflows

---

### ❌ **2. No Variable Name Storage**

- Saved data has no associated variable name
- User must manually assign loaded data to a variable
- Cannot reconstruct workspace like MATLAB's `load()`

**MATLAB behavior**:
```matlab
save('test.mat', 'myData')
clear
load('test.mat')  % Variable 'myData' appears in workspace
```

**Current behavior**:
```matlab
SaveMatFast('test.bin', myData)
clear
data = LoadMatFast('test.bin')  % Must manually name it
```

---

### ⚠️ **3. Fixed Comment Size** (Line 21-28)

```matlab
if numel(comment) > 42
    comment = comment(1:42);  % Truncation
end
```

- Comment limited to exactly 42 characters
- Padding with spaces is wasteful
- Should use variable-length field

---

### ⚠️ **4. No Format Version**

- No version number in header
- Future format changes will break compatibility
- No way to detect old vs new format

---

### ⚠️ **5. Limited Error Handling**

**Missing validations**:
- No check for write permissions before creating file
- No verification of free disk space for large files
- Minimal error messages (e.g., line 16: "Failed to open file")
- No recovery mechanism for partial writes

---

### ⚠️ **6. Header Size Calculation** (Line 74)

```matlab
header.headersize = 2 + 1 + 14 + 8 + 1 + 6 + 8 + header.dataNumberDimensions*8 + 1 + 42;
```

**Issues**:
- Magic numbers without explanation
- Brittle - easy to break if header changes
- Should be calculated from actual written data

---

## Code Quality Issues

### **SaveMatFast.m**

#### Line 77: Typo in error message
```matlab
error(['Error openning file: ', errormsg])  % "openning" -> "opening"
```

#### Line 21-28: Awkward comment padding
```matlab
tmp = repmat(' ',1,42);
tmp(1:numel(comment)) = comment;
comment = tmp;
```
Should use `pad()` or sprintf formatting

#### Line 109-113: Redundant type checking
```matlab
if strcmp(datatype, 'logical')
    fwrite(fp, uint8(datamat), 'uint8');
else
    fwrite(fp, datamat, datatype);
end
```
Could be handled in single path

### **LoadMatFast.m**

#### Line 11: Weak error message
```matlab
error([filepath,' does not exist.'])
```
Should suggest checking path or listing available files

#### Line 53-56: Checksum validation
```matlab
if checksum(1)~=1 || checksum(2) ~= 1
    fclose(fp);
    error('Check sum mismatch!')
```
- Inconsistent spacing around operators
- Trivial checksum (always [1,1]) - not useful
- Should use CRC32 or MD5 for real verification

#### Line 138-161: Inefficient struct loading
Pre-allocates cell array then converts to struct. Could be optimized.

---

## Performance Analysis

### **Advantages over MATLAB save/load**

1. **No compression overhead** - Direct binary write
2. **No variable name lookup** - Single variable is faster
3. **Optimized for large matrices** - Sequential write

### **Disadvantages**

1. **No compression** - Large files for sparse/redundant data
2. **No lazy loading** - Entire file must be read
3. **No partial access** - Can't load subset of data
4. **Recursive overhead** - Structs/cells have per-element overhead

---

## Security Concerns

### ⚠️ **File Overwrite** (Line 76)
```matlab
[fp,errormsg]= fopen(filepath,'w');  % Overwrites without warning
```
Should check if file exists and prompt user

### ⚠️ **No Validation of Loaded Data**
- File could be corrupted or malicious
- No bounds checking on dimensions
- Could cause memory exhaustion with crafted file

### ✅ **Type Safety**
- Type IDs prevent type confusion
- Dimension validation exists

---

## Compatibility Issues

### **Endianness**
- No explicit endianness handling
- Will fail when transferring files between big-endian and little-endian systems
- MATLAB save handles this automatically

### **Path Separators**
- Uses platform-specific file paths
- No cross-platform path handling

---

## Recommendations

### **Priority 1: Critical**

1. ✅ **Add multi-variable support**
   - Store variable names in header
   - Support syntax: `SaveMatFast('file.bin', 'var1', var1, 'var2', var2)`
   - Load returns struct with field names matching variable names

2. ✅ **Add format version**
   - Include version number in header
   - Allow backwards compatibility checking

3. ✅ **Improve error handling**
   - Validate inputs before writing
   - Check disk space availability
   - Better error messages with context

### **Priority 2: High**

4. **Fix comment field**
   - Use variable-length comment
   - Remove 42-character limit

5. **Add real checksum**
   - Replace [1,1] with CRC32 or MD5
   - Validate file integrity on load

6. **Add compression option**
   - Optional zlib/gzip compression
   - Trade speed for file size

### **Priority 3: Medium**

7. **Add partial loading**
   - Load specific variables by name
   - Memory map for very large files

8. **Performance optimization**
   - Batch write operations
   - Optimize struct serialization

9. **Cross-platform compatibility**
   - Handle endianness
   - Platform-independent paths

---

## Test Coverage (Currently: 0%)

### **Missing Tests**

- ❌ No unit tests
- ❌ No integration tests
- ❌ No performance benchmarks
- ❌ No edge case validation
- ❌ No corruption recovery tests

### **Recommended Test Suite**

1. **Basic I/O Tests**
   - Save/load each supported type
   - Round-trip verification
   - Size verification

2. **Complex Data Tests**
   - Nested structs (5+ levels deep)
   - Large cell arrays (1M+ elements)
   - Mixed-type structures

3. **Edge Cases**
   - Empty arrays
   - Scalar values
   - Multi-dimensional arrays (3D, 4D, etc.)
   - Special values (NaN, Inf, -Inf)

4. **Error Handling Tests**
   - Invalid file paths
   - Corrupted files
   - Insufficient disk space
   - Permission errors

5. **Performance Tests**
   - Benchmark vs MATLAB save/load
   - Large file handling (>1GB)
   - Memory usage profiling

---

## Conclusion

**Current State**: Functional but incomplete implementation suitable for single-variable, simple use cases.

**Recommended Action**: Refactor to support multiple variables and add comprehensive testing before production use.

**Risk Assessment**:
- **LOW** for single-variable numeric arrays
- **MEDIUM** for complex nested structures
- **HIGH** for production use without multi-variable support

---

## Refactoring Plan

See refactored versions:
- `SaveFast.m` - Multi-variable support, better error handling
- `LoadFast.m` - Multi-variable support, workspace integration
- `test_MatFast.m` - Comprehensive test suite
