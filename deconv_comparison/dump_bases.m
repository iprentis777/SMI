function dump_bases()
% dump_bases()
%
% Write SMI's even spherical harmonic basis on the shared evaluation
% directions, with and without the Condon-Shortley phase, so that
% check_conventions.py can compare both against MRtrix's basis. Run after
% setup_protocol.py.
more off
IO = binio();
ev = IO.load('eval_dirs');
IO.save('Y_smi_cs0', SMI.get_even_SH(ev, 6, 0));
IO.save('Y_smi_cs1', SMI.get_even_SH(ev, 6, 1));
fprintf('wrote Y_smi_cs0 and Y_smi_cs1 [%d x %d]\n', size(ev,1), 28);
end
