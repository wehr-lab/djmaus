%ProcessSession - backward-compatible SCRIPT wrapper around the function ProcessSession2.
%Run from a Bonsai folder. Calls ProcessSession2(pwd) and copies every output into the
%caller's workspace under the same variable names the old script created (BonsaiPath,
%EphysPath, EphysPath_KS, SortedUnitsFile, stimtrace, timestamps, ...), so existing code
%that does "ProcessSession" and then uses those variables keeps working. All outputs are
%(re)assigned on every call, so nothing is carried over from a previous session.
%New code should call PS = ProcessSession2(BonsaiPath, ...) directly.
%The original script is preserved as ProcessSession_legacy_script.m.
%-mike/claude 10.2.26
ProcessSession_PS__=ProcessSession2(pwd);
for ProcessSession_f__=fieldnames(ProcessSession_PS__)'
    eval([ProcessSession_f__{1} '=ProcessSession_PS__.(ProcessSession_f__{1});']);
end
clear ProcessSession_PS__ ProcessSession_f__
