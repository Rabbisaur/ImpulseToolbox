function [userview systemview] = GetMemory
% Emulate the behavior of function memory on windows system
% tested on Ubuntu 11.04
% Created: Jan.4, 2012  Wang Feng@BNU
if isunix
    %% get MATLAB memory usage values
    if ~ismac
        [~,page_size] = unix('getconf PAGE_SIZE');
        ratio = str2double(page_size) / 1024;
        [~, a] = unix('ps a|grep MATLAB');
        matlabpid = textscan(a,' %d',1);
        cmd = ['cat /proc/' num2str(cell2mat(matlabpid)) '/statm'];
        [~, a] = unix(cmd);
        b = textscan(a,'%f');
        if numel(b{1})<2
            matlab_mem = -1;
        else
            matlab_mem = double(b{1}(2))*ratio*1024;
        end
    else
        matlab_mem = -1;
    end
    %% get vmstat values
    if ~ismac
        cmd = 'cat /proc/meminfo';
        [~, b] = unix(cmd);
        a = textscan(b,'%s');
        a = a{1};
        idx = strcmpi(a,'MemTotal:');
        idx = find(idx);
        total_memory = str2double(a(idx+1))*1024;
%         free_memory = str2double(a{1}(5))*1024;
%         buffers_memory = str2double(a{1}(8))*1024;
%         cached_memory =str2double(a{1}(11))*1024;
%         total_swap = str2double(a{1}(41))*1024;
%         free_swap = str2double(a{1}(44))*1024;
        % availiable physical memory = free + buffer + cached
%         available_memory = free_memory + buffers_memory + cached_memory;
        idx = strcmpi(a,'MemAvailable:');
        idx = find(idx);
        available_memory = str2double(a(idx+1))*1024;
        
        idx = strcmpi(a,'SwapFree:');
        idx = find(idx);
        free_swap = str2double(a(idx+1))*1024;
        
        idx = strcmpi(a,'SwapTotal:');
        idx = find(idx);
        total_swap = str2double(a(idx+1))*1024;
    else
        % % currently no implementation
        error('Sorry, funtion GetMemory currently does not know know to work on Mac OS')
    end
    %% get userview values
    % Max Possible Array Bytes
    userview.MaxPossibleArrayBytes = available_memory + free_swap;
    %MemAvailable All Arrays
    userview.MemAvailableAllArrays = available_memory + free_swap;
    %MemUsed MATLAB
    userview.MemUsedMATLAB = matlab_mem;
    %% get systemview values
    % Virtual Address Space
    systemview.VirtualAddressSpace.Available = free_swap;
    systemview.VirtualAddressSpace.Total = total_swap;
    % systemview
    %System Memory
    systemview.SystemMemory.Available = available_memory + free_swap;
    %Physical Memory
    systemview.PhysicalMemory.Available = available_memory;
    systemview.PhysicalMemory.Total = total_memory;
elseif ispc
    [userview systemview] = memory;
end