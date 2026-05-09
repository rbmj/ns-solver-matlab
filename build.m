cfg = coder.gpuConfig('mex');
cfg.GenerateReport = true;
cfg.GpuConfig.CompilerFlags = '--fmad=false --compiler-bindir=/usr/local/MATLAB/R2025b/TOOLCHAIN';

codegen -config cfg recalculate_parameters
codegen -config cfg advance_time