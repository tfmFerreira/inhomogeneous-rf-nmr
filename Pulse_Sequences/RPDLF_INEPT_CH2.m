function [] = WriteSimpsonIn_RPDLF(varargin)
%WriteSimpsonIn_RPDLF: Writes an .in file for an RPDLF experiment to be interpreted by SIMPSON.


l1=9;

% evaluate user input
default_np = 0;
default_dipole = [1 3 0 0 0 0];
default_dipoleH = [2 3 0 0 0 0];
default_dipoleHH = [1 2 0 0 0 0];
default_jcoupling = [1 3 140 0 0 0 0 0];
default_jcouplingb = [2 3 140 0 0 0 0 0];
default_acc = 0;
default_accb = 0;
default_dip1 = 0;
default_dip2 = 0;
default_dipHH = 0;
default_shiftH = 0;
default_shiftC = 0;
default_CSAH = 0;
default_CSAC=0;
default_phasecount=0;
default_rep=0;
default_gamma=0;
default_MAS=0;
default_maxdt=0;

p= inputParser;

np='NP';
dipole='dipole';
dipoleH='dipoleH';
dipoleHH='dipoleHH';
jcoupling='jcoupling';
jcouplingb='jcoupling';
acc='acc';
accb='accb';
dip1='dip1';
dip2='dip2';
dipHH='dipHH';
phasecount='phasecount';
shiftH='shiftH';
shiftC='shiftC';
CSAH='CSAH';
CSAC='CSAC';
rep='rep';
gamma='gamma';
MAS='MAS';
maxdt='maxdt';

addParameter(p,'NP',default_np);
addParameter(p,'dipole',default_dipole);
addParameter(p,'dipoleH',default_dipoleH);
addParameter(p,'dipoleHH',default_dipoleHH);
addParameter(p,'jcoupling',default_jcoupling)
addParameter(p,'jcouplingb',default_jcouplingb)
addParameter(p,'acc',default_acc)
addParameter(p,'accb',default_accb)
addParameter(p,'dip1',default_dip1)
addParameter(p,'dip2',default_dip2)
addParameter(p,'dipHH',default_dipHH)
addParameter(p,'phasecount',default_phasecount)
addParameter(p,'shiftH',default_shiftH)
addParameter(p,'shiftC',default_shiftC)
addParameter(p,'CSAH',default_CSAH)
addParameter(p,'CSAC',default_CSAC)
addParameter(p,'rep',default_rep)
addParameter(p,'gamma',default_gamma)
addParameter(p,'MAS',default_MAS)
addParameter(p,'maxdt',default_maxdt)

parse(p,varargin{:});

np = p.Results.NP;
dipole = p.Results.dipole;
dipoleH = p.Results.dipoleH;
dipoleHH = p.Results.dipoleHH;
jcoupling = p.Results.jcoupling;
jcouplingb = p.Results.jcouplingb;
acc = p.Results.acc;
acc = 1.0-acc*0.01;
accb = p.Results.accb;
accb = 1.0-accb*0.01;
dip1=p.Results.dip1;
dip2=p.Results.dip2;
dipHH=p.Results.dipHH;
phasecount=p.Results.phasecount;
shiftH=p.Results.shiftH;
shiftC=p.Results.shiftC;
CSAH=p.Results.CSAH;
CSAC=p.Results.CSAC;
rep=p.Results.rep;
gamma=p.Results.gamma;
MAS=p.Results.MAS;
maxdt=p.Results.maxdt;

ph1=[1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3]+1;
ph2=[1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3]+1;
ph3=[0 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2 0 0 2 2]+1;
ph4=ph2;
ph5=[1 1 1 1 2 2 2 2 3 3 3 3 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 0 0 0 0]+1;
ph6=[1 3 1 3 2 0 2 0 1 3 1 3 2 0 2 0 1 3 1 3 2 0 2 0 1 3 1 3 2 0 2 0]+1;
ph7=[70.0];
ph8=[290.0];
ph9=[250.0];
ph10=[110.0];
ph11=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]+1;
ph31=[3 3 1 1 0 0 2 2 3 3 1 1 0 0 2 2 3 3 1 1 0 0 2 2 3 3 1 1 0 0 2 2]+1;

p_phase=[0.0 90.0 180.0 270.0];
d_phase={'I3x' 'I3y' '-I3x' '-I3y'};


% ----------write SIMPSON .in file----------------

% create new (temporary) file

file=strcat('../TMP/tmp.in');
fID=fopen(file,'w');
fprintf(fID,'%s\n','# (R18)^1_7');
fprintf(fID,'%s\n',' ');

% spinsystem
fprintf(fID,'%s\n','spinsys {');
fprintf(fID,'%4s%6s%3s%3s%4s\n',' ','nuclei','1H','1H','13C');
fprintf(fID,'%4s%7s%3s%4s\n',' ','channels','1H','13C');
fprintf(fID,'%4s%6s%2s%2s%8s%4s%4s%4s\n',' ','dipole',mat2str(dipole(1)),mat2str(dipole(2)),mat2str(dip1),mat2str(dipole(4)),mat2str(dipole(5)),mat2str(dipole(6)));
fprintf(fID,'%4s%6s%2s%2s%8s%4s%4s%4s\n',' ','dipole',mat2str(dipoleH(1)),mat2str(dipoleH(2)),mat2str(dip2),mat2str(dipoleH(4)),mat2str(dipoleH(5)),mat2str(dipoleH(6)));
fprintf(fID,'%4s%6s%2s%2s%8s%4s%4s%4s\n',' ','dipole',mat2str(dipoleHH(1)),mat2str(dipoleHH(2)),mat2str(dipHH),mat2str(dipoleHH(4)),mat2str(dipoleHH(5)),mat2str(dipoleHH(6)));
fprintf(fID,'%4s%6s%6s%1s%3.2f%2s%3.2f%1s%4s%4s%4s%4s\n',' ','shift',mat2str(dipoleHH(1)),' ',shiftH,'p ',CSAH,'p','0.5','0','0','0');
fprintf(fID,'%4s%6s%6s%1s%3.2f%2s%3.2f%1s%4s%4s%4s%4s\n',' ','shift',mat2str(dipoleHH(2)),' ',shiftH,'p ',CSAH,'p','0.5','0','0','0');
fprintf(fID,'%4s%6s%6s%1s%3.2f%2s%3.2f%1s%4s%4s%4s%4s\n',' ','shift',mat2str(dipole(2)),' ',shiftC,'p ',CSAC,'p','0.5','0','0','0');
fprintf(fID,'%4s%6s%2s%2s%5s%2s%2s%2s%2s%2s\n',' ','jcoupling',mat2str(jcoupling(1)),mat2str(jcoupling(2)),mat2str(jcoupling(3)),mat2str(jcoupling(4)),mat2str(jcoupling(5)),mat2str(jcoupling(6)),mat2str(jcoupling(7)),mat2str(jcoupling(8)));
fprintf(fID,'%4s%6s%2s%2s%5s%2s%2s%2s%2s%2s\n',' ','jcoupling',mat2str(jcouplingb(1)),mat2str(jcouplingb(2)),mat2str(jcoupling(3)),mat2str(jcoupling(4)),mat2str(jcoupling(5)),mat2str(jcoupling(6)),mat2str(jcoupling(7)),mat2str(jcoupling(8)));
fprintf(fID,'%1s\n','}');

% parameters
fprintf(fID,'%s\n','par {');
fprintf(fID,'%4s%16s%10s\n',' ','proton_frequency','400e6');        

if rep > 99
fprintf(fID,'%4s%12s%4s%3s\n',' ','crystal_file','rep',mat2str(rep));
elseif rep >999
fprintf(fID,'%4s%12s%4s%4s\n',' ','crystal_file','rep',mat2str(rep));
else
end

fprintf(fID,'%4s%12s%14s\n',' ','dipole_check','false');
fprintf(fID,'%4s%12s%10s\n',' ','gamma_angles',num2str(gamma));
fprintf(fID,'%4s%9s%16s\n',' ','spin_rate',num2str(MAS));
fprintf(fID,'%4s%11s%11s\n',' ','variable lb','0');
fprintf(fID,'%4s%2s%32s\n',' ','sw','spin_rate/2.0');
fprintf(fID,'%4s%7s%18s\n',' ','verbose','1101');
fprintf(fID,'%4s%13s%18s\n',' ','conjugate_fid','true');

fprintf(fID,'%4s%2s%22s\n',' ','np',mat2str(np));                     % number of points in direct dimension
fprintf(fID,'%4s%14s%10s\n',' ','start_operator','Inz');
fprintf(fID,'%4s%15s%10s\n',' ','detect_operator',d_phase{ph31(phasecount)});
fprintf(fID,'%s\n',' ');
fprintf(fID,'%1s\n','}');
fprintf(fID,'%s\n',' ');
fprintf(fID,'%s\n',' ');

% pulse sequence
% parameters
fprintf(fID,'%s\n','proc pulseq {} {');
fprintf(fID,'%8s%s\n',' ','global par');
fprintf(fID,'%1s\n',' ');
fprintf(fID,'%8s%4.2f\n','maxdt ',maxdt);     % max timestep over which H is time-independent
fprintf(fID,'%1s\n',' ');
fprintf(fID,'%8s%s%3.1f%s\n',' ','set rf2 [expr ',l1,'*$par(spin_rate)]');
fprintf(fID,'%8s%s\n',' ','set rf5 78125.0');
fprintf(fID,'%8s%s%4.2f%s\n',' ','set rf4 [expr ',acc,'*$rf5]'); % pulse with accuracey 'acc'
fprintf(fID,'%8s%s%4.2f%s\n',' ','set rf3 [expr ',accb,'*$rf5]'); % pulse with accuracey 'acc'
fprintf(fID,'%8s%s%4.2f%1s%3.1f%s\n',' ','set rf [expr ',acc,'*',l1,'*$par(spin_rate)]'); % pulse with accuracey 'acc'
fprintf(fID,'%8s%s%4.2f%1s%3.1f%s\n',' ','set rf6 [expr ',accb,'*',2*l1,'*$par(spin_rate)]'); % pulse with accuracey 'acc'
fprintf(fID,'%8s%s\n',' ','set t180 [expr 0.5e6/$rf2]');
fprintf(fID,'%8s%s\n',' ','set t90 [expr 0.25e6/$rf2]');
fprintf(fID,'%8s%s\n',' ','set t45 [expr 0.125e6/$rf2]');
fprintf(fID,'%8s%s\n',' ','set t180I [expr 0.5e6/$rf5]');
fprintf(fID,'%8s%s\n',' ','set t90I [expr 0.25e6/$rf5]');
fprintf(fID,'%8s%7s%5.2f\n',' ','set ph1 ',p_phase(ph1(phasecount)));
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph2 ',p_phase(ph2(phasecount)));
fprintf(fID,'%8s%7s%5.2f\n',' ','set ph3 ',p_phase(ph3(phasecount)));
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph4 ',p_phase(ph4(phasecount)));
fprintf(fID,'%8s%7s%5.2f\n',' ','set ph5 ',p_phase(ph5(phasecount)));
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph6 ',p_phase(ph6(phasecount)));
fprintf(fID,'%8s%7s%5.2f\n',' ','set ph7 ',ph7);
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph8 ',ph8);
fprintf(fID,'%8s%7s%5.2f\n',' ','set ph9 ',ph9);
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph10 ',ph10);
fprintf(fID,'%8s%8s%5.2f\n',' ','set ph11 ',p_phase(ph11(phasecount)));
fprintf(fID,'%8s%s\n',' ','set d1 [expr 0.0012e6]');
fprintf(fID,'%8s%s\n',' ','set d2 [expr 0.002e6]');

% reset current propagator
fprintf(fID,'%s\n','reset');
fprintf(fID,'%5s\n',' ');


% pulses

% first R-block
for i=1:l1
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph7 0 $ph2');
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph8 0 $ph2');
end
fprintf(fID,'%8s%s\n',' ','store 1'); % store propagator
fprintf(fID,'%s\n','reset');
fprintf(fID,'%s\n',' ');
% 2nd R-block
for i=1:l1
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph9 0 $ph2');
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph10 0 $ph2');
end
fprintf(fID,'%8s%s\n',' ','store 2'); % store propagator
fprintf(fID,'%5s\n',' ');
fprintf(fID,'%s\n',' ');

% refocused INEPT
fprintf(fID,'%s\n','reset');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph1 0 y'); % ideal pulse
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph4');
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph3 $rf3 $ph5');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph6');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','acq');
fprintf(fID,'%s\n',' ');
fprintf(fID,'%s\n',' ');

% middle of R-block sequence: R-180-R
fprintf(fID,'%s\n','reset');
for i=1:l1-1
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph7 0 $ph2');
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph8 0 $ph2');
end
fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph7 0 $ph2');
fprintf(fID,'%8s%s\n',' ','pulse $t90 $rf $ph8 0 $ph2');
fprintf(fID,'%8s%s\n',' ','pulse $t45 $rf $ph8 0 $ph11');
fprintf(fID,'%8s%s\n',' ','pulse $t45 $rf $ph8 $rf6 $ph11');
fprintf(fID,'%8s%s\n',' ','pulse $t45 $rf $ph9 $rf6 $ph11');
fprintf(fID,'%8s%s\n',' ','pulse $t45 $rf $ph9 0 $ph11');
fprintf(fID,'%8s%s\n',' ','pulse $t90 $rf $ph9 0 $ph2');
fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph10 0 $ph2');
for i=1:l1-1
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph9 0 $ph2');
    fprintf(fID,'%8s%s\n',' ','pulse $t180 $rf $ph10 0 $ph2');
end
fprintf(fID,'%s\n','store 3');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph1 0 y'); % ideal pulse
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph4');
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph3 $rf3 $ph5');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph6');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','acq');
fprintf(fID,'%s\n',' ');
fprintf(fID,'%s\n',' ');

%
fprintf(fID,'%8s%s\n',' ','for {set j 2} {$j < $par(np)} {incr j} {');
fprintf(fID,'%8s%s\n',' ','reset');
fprintf(fID,'%8s%s\n',' ','prop 1'); % first R-blocks
fprintf(fID,'%8s%s\n',' ','prop 3'); % R-180-R
fprintf(fID,'%8s%s\n',' ','prop 2'); % second R-blocks
fprintf(fID,'%8s%s\n',' ','store 3');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph1 0 y'); % ideal pulse
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph4');
fprintf(fID,'%s\n','delay $d2');
fprintf(fID,'%s\n','pulse $t90I $rf4 $ph3 $rf3 $ph5');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','pulse $t180I $rf4 $ph2 $rf3 $ph6');
fprintf(fID,'%s\n','delay $d1');
fprintf(fID,'%s\n','acq}');
fprintf(fID,'%s\n','}');
fprintf(fID,'%s\n',' ');
fprintf(fID,'%s\n',' ');

% processing section

fprintf(fID,'%s\n','proc main {} {');
fprintf(fID,'%4s%s\n',' ','global par');
fprintf(fID,'%4s%s\n',' ','set f [fsimpson]');
fprintf(fID,'%4s%s\n',' ','fsave $f $par(name).fid');
fprintf(fID,'%4s%s\n',' ','fsave $f $par(name)_TD.dat -xreim'); % direct dim
fprintf(fID,'%4s%s\n',' ','fft $f');
fprintf(fID,'%4s%s\n',' ','fsave $f $par(name).spe');
fprintf(fID,'%4s%s\n',' ','fsave $f $par(name)_FT.dat -xreim');
fprintf(fID,'%s\n','}');

fclose(fID);


end



