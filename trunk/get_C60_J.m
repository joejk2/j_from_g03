disp ("Assumes the presence of the following files: MOs_1.txt, MOs_2.txt, MOs_pair.txt, S_1.txt, S_2.txt, S_pair.txt, Evls_pair.txt")
disp ("!!! Octave counts from 1 !!!")

%nbasis   =input("Number of MOs (of monomer)");		%... of monomer
%nhomo_mon=input("HOMO (of monomer)");		%... of monomer (!! Octave counts from 1 !!)
%nbasis=446
%nhomo_mon=73
%Naphthalene
%nbasis=204
%nhomo_mon=34

%C60
nbasis=900
nhomo_mon=180

%sprintf('read data')
load('MOs_1.txt');
load('MOs_2.txt');
load('MOs_pair.txt');
load('S_1.txt');
load('S_2.txt');
load('S_pair.txt');
load('Evls_pair.txt');

zero_mat=zeros(nbasis);
MOs_mons = [MOs_1, zero_mat; zero_mat, MOs_2];
S_mons   = [S_1,   zeros(nbasis); zeros(nbasis), S_2 ];
D_mons   = sqrtm(S_mons);
D_pair   = sqrtm(S_pair);

%b3lyp - have to orthogonalise
MOs_mons_orth = D_mons * MOs_mons;
MOs_pair_orth = D_pair * MOs_pair;

B = MOs_mons_orth' * MOs_pair_orth;
Evls=diag(Evls_pair');
H_eff = B * Evls * B';

%test
%B_test = MOs_mons' * MOs_pair;
%H_eff_test = B_test * Evls * B_test';
%H00_test=H_eff_test(nhomo_mon+nbasis, nhomo_mon)
%           -0.5*(H_eff_test(nhomo_mon+nbasis,nhomo_mon+nbasis)+H_eff_test(nhomo_mon,nhomo_mon))*S_pair(nhomo_mon+nbasis,nhomo_mon)
%	  )/(1-S_pair(nhomo_mon+nbasis,nhomo_mon)*S_pair(nhomo_mon+nbasis,nhomo_mon))

%<HOMO | F | HOMO > 
%<LUMO | F | LUMO > 
%Output.  N.B. Energies already converted to eV by rewrite_S_phi_E.cpp

%non-degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%E1_E2=(H_eff(nhomo_mon,nhomo_mon)-H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis))
%half_HOMO_diff=(Evls(2*nhomo_mon,2*nhomo_mon)-Evls(2*nhomo_mon-1,2*nhomo_mon-1))/2
%Non_degenerate_HOMO_coupling=H_eff(nhomo_mon+nbasis,   nhomo_mon  )
%Non_degenerate_LUMO_coupling=H_eff(nhomo_mon+1+nbasis, nhomo_mon+1)
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%Non_degenerate_LUMO_coupling=abs(L00)

%E HOMO A
E_H0_A=H_eff(nhomo_mon,   nhomo_mon  )
E_H1_A=H_eff(nhomo_mon-1, nhomo_mon-1)
E_H2_A=H_eff(nhomo_mon-2, nhomo_mon-2)
E_H3_A=H_eff(nhomo_mon-3, nhomo_mon-3)
E_H4_A=H_eff(nhomo_mon-4, nhomo_mon-4)
%E HOMO B
E_H0_B=H_eff(nhomo_mon+nbasis,   nhomo_mon+nbasis  )
E_H1_B=H_eff(nhomo_mon+nbasis-1, nhomo_mon+nbasis-1)
E_H2_B=H_eff(nhomo_mon+nbasis-2, nhomo_mon+nbasis-2)
E_H3_B=H_eff(nhomo_mon+nbasis-3, nhomo_mon+nbasis-3)
E_H4_B=H_eff(nhomo_mon+nbasis-4, nhomo_mon+nbasis-4)

%Difference in HOMO energies
avg_E_HOMO_A_minus_E_HOMO_B=(E_H0_A + E_H1_A + E_H2_A + E_H3_A + E_H4_A - E_H0_B - E_H1_B - E_H2_B - E_H3_B - E_H4_B)/5

%E LUMO A
E_L0_A=H_eff(nhomo_mon+1, nhomo_mon+1)
E_L1_A=H_eff(nhomo_mon+2, nhomo_mon+2)
E_L2_A=H_eff(nhomo_mon+3, nhomo_mon+3)
E_L3_A=H_eff(nhomo_mon+4, nhomo_mon+4)
E_L4_A=H_eff(nhomo_mon+5, nhomo_mon+5)
%E LUMO B
E_L0_B=H_eff(nhomo_mon+nbasis+1, nhomo_mon+nbasis+1)
E_L1_B=H_eff(nhomo_mon+nbasis+2, nhomo_mon+nbasis+2)
E_L2_B=H_eff(nhomo_mon+nbasis+3, nhomo_mon+nbasis+3)
E_L3_B=H_eff(nhomo_mon+nbasis+4, nhomo_mon+nbasis+4)
E_L4_B=H_eff(nhomo_mon+nbasis+5, nhomo_mon+nbasis+5)

%Difference in LUMO energies
avg_E_LUMO_A_minus_E_LUMO_B=(E_L0_A+E_L1_A+E_L2_A+E_L3_A+E_L4_A-E_L0_B-E_L1_B-E_L2_B-E_L3_B-E_L4_B)/5

%tests
should_be_zero=H_eff(nhomo_mon, nhomo_mon-1)
should_be_zero=H_eff(nhomo_mon, nhomo_mon-2)
should_be_zero=H_eff(nhomo_mon, nhomo_mon-3)
should_be_zero=H_eff(nhomo_mon, nhomo_mon-4)
should_be_zero=H_eff(nhomo_mon, nhomo_mon-5)
should_be_zero=H_eff(nhomo_mon+nbasis+4, nhomo_mon+nbasis+0)
should_be_zero=H_eff(nhomo_mon+nbasis+3, nhomo_mon+nbasis+1)
should_be_zero=H_eff(nhomo_mon+nbasis+2, nhomo_mon+nbasis+2)
should_be_zero=H_eff(nhomo_mon+nbasis+1, nhomo_mon+nbasis+3)
should_be_zero=H_eff(nhomo_mon+nbasis+0, nhomo_mon+nbasis+4)

%H0
H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  )
H01=H_eff(nhomo_mon+nbasis,   nhomo_mon-1)
H02=H_eff(nhomo_mon+nbasis,   nhomo_mon-2)
H03=H_eff(nhomo_mon+nbasis,   nhomo_mon-3)
H04=H_eff(nhomo_mon+nbasis,   nhomo_mon-4)
%H1
H10=H_eff(nhomo_mon+nbasis-1, nhomo_mon  )
H11=H_eff(nhomo_mon+nbasis-1, nhomo_mon-1)
H12=H_eff(nhomo_mon+nbasis-1, nhomo_mon-2)
H13=H_eff(nhomo_mon+nbasis-1, nhomo_mon-3)
H14=H_eff(nhomo_mon+nbasis-1, nhomo_mon-4)
%H2
H20=H_eff(nhomo_mon+nbasis-2, nhomo_mon  )
H21=H_eff(nhomo_mon+nbasis-2, nhomo_mon-1)
H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-2)
H23=H_eff(nhomo_mon+nbasis-2, nhomo_mon-3)
H24=H_eff(nhomo_mon+nbasis-2, nhomo_mon-4)
%H3
H30=H_eff(nhomo_mon+nbasis-3, nhomo_mon  )
H31=H_eff(nhomo_mon+nbasis-3, nhomo_mon-1)
H32=H_eff(nhomo_mon+nbasis-3, nhomo_mon-2)
H33=H_eff(nhomo_mon+nbasis-3, nhomo_mon-3)
H34=H_eff(nhomo_mon+nbasis-3, nhomo_mon-4)
%H4
H40=H_eff(nhomo_mon+nbasis-4, nhomo_mon  )
H41=H_eff(nhomo_mon+nbasis-4, nhomo_mon-1)
H42=H_eff(nhomo_mon+nbasis-4, nhomo_mon-2)
H43=H_eff(nhomo_mon+nbasis-4, nhomo_mon-3)
H44=H_eff(nhomo_mon+nbasis-4, nhomo_mon-4)


%HOMO couplings
avg_5fold_degenerate_HOMO_coupling=( H00 + H01 + H02 + H03 + H04
				   + H10 + H11 + H12 + H13 + H14
				   + H20 + H21 + H22 + H23 + H24
				   + H30 + H31 + H32 + H33 + H34
				   + H40 + H41 + H42 + H43 + H44 )/25
rms_5fold_degenerate_HOMO_coupling=sqrt( H00*H00 + H01*H01 + H02*H02 + H03*H03 + H04*H04
				       + H10*H10 + H11*H11 + H12*H12 + H13*H13 + H14*H14
				       + H20*H20 + H21*H21 + H22*H22 + H23*H23 + H24*H24
				       + H30*H30 + H31*H31 + H32*H32 + H33*H33 + H34*H34
				       + H40*H40 + H41*H41 + H42*H42 + H43*H43 + H44*H44 )/25


%L0
L00=H_eff(nhomo_mon+nbasis+1,   nhomo_mon+1)
L01=H_eff(nhomo_mon+nbasis+1,   nhomo_mon+2)
L02=H_eff(nhomo_mon+nbasis+1,   nhomo_mon+3)
L03=H_eff(nhomo_mon+nbasis+1,   nhomo_mon+4)
L04=H_eff(nhomo_mon+nbasis+1,   nhomo_mon+5)
%L1
L10=H_eff(nhomo_mon+nbasis+2, nhomo_mon+1)
L11=H_eff(nhomo_mon+nbasis+2, nhomo_mon+2)
L12=H_eff(nhomo_mon+nbasis+2, nhomo_mon+3)
L13=H_eff(nhomo_mon+nbasis+2, nhomo_mon+4)
L14=H_eff(nhomo_mon+nbasis+2, nhomo_mon+5)
%L2
L20=H_eff(nhomo_mon+nbasis+3, nhomo_mon+1)
L21=H_eff(nhomo_mon+nbasis+3, nhomo_mon+2)
L22=H_eff(nhomo_mon+nbasis+3, nhomo_mon+3)
L23=H_eff(nhomo_mon+nbasis+3, nhomo_mon+4)
L24=H_eff(nhomo_mon+nbasis+3, nhomo_mon+5)
%L3
L30=H_eff(nhomo_mon+nbasis+4, nhomo_mon+1)
L31=H_eff(nhomo_mon+nbasis+4, nhomo_mon+2)
L32=H_eff(nhomo_mon+nbasis+4, nhomo_mon+3)
L33=H_eff(nhomo_mon+nbasis+4, nhomo_mon+4)
L34=H_eff(nhomo_mon+nbasis+4, nhomo_mon+5)
%L4
L40=H_eff(nhomo_mon+nbasis+5, nhomo_mon+1)
L41=H_eff(nhomo_mon+nbasis+5, nhomo_mon+2)
L42=H_eff(nhomo_mon+nbasis+5, nhomo_mon+3)
L43=H_eff(nhomo_mon+nbasis+5, nhomo_mon+4)
L44=H_eff(nhomo_mon+nbasis+5, nhomo_mon+5)

%LUMO couplings
avg_5fold_degenerate_LUMO_coupling=( L00 + L01 + L02 + L03 + L04
				   + L10 + L11 + L12 + L13 + L14
				   + L20 + L21 + L22 + L23 + L24
				   + L30 + L31 + L32 + L33 + L34
				   + L40 + L41 + L42 + L43 + L44 )/25
rms_5fold_degenerate_LUMO_coupling=sqrt( L00*L00 + L01*L01 + L02*L02 + L03*L03 + L04*L04
				       + L10*L10 + L11*L11 + L12*L12 + L13*L13 + L14*L14
				       + L20*L20 + L21*L21 + L22*L22 + L23*L23 + L24*L24
				       + L30*L30 + L31*L31 + L32*L32 + L33*L33 + L34*L34
				       + L40*L40 + L41*L41 + L42*L42 + L43*L43 + L44*L44 )/25

%doubly degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%H01=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%H10=H_eff(nhomo_mon+nbasis-1, nhomo_mon  );
%H11=H_eff(nhomo_mon+nbasis-1, nhomo_mon-1);
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L01=H_eff(nhomo_mon+nbasis+1, nhomo_mon);
%L10=H_eff(nhomo_mon+nbasis,   nhomo_mon+1);
%L11=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%Doubly_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H10)+abs(H11))/4
%Doubly_degenerate_LUMO_coupling=(abs(L00)+abs(L01)+abs(L10)+abs(L11))/4

%triply degenerate
%H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%H01=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%H02=H_eff(nhomo_mon+nbasis,   nhomo_mon-2);
%H10=H_eff(nhomo_mon+nbasis-1, nhomo_mon  );
%H11=H_eff(nhomo_mon+nbasis-1, nhomo_mon-1);
%H12=H_eff(nhomo_mon+nbasis-1, nhomo_mon-2);
%H20=H_eff(nhomo_mon+nbasis-2, nhomo_mon  );
%H21=H_eff(nhomo_mon+nbasis-2, nhomo_mon-1);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-2);
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L01=H_eff(nhomo_mon+nbasis+1, nhomo_mon);
%L02=H_eff(nhomo_mon+nbasis+1, nhomo_mon-1);
%L10=H_eff(nhomo_mon+nbasis,   nhomo_mon+1);
%L11=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
%L12=H_eff(nhomo_mon+nbasis,   nhomo_mon-1);
%L20=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%L21=H_eff(nhomo_mon+nbasis+1, nhomo_mon  );
%L22=H_eff(nhomo_mon+nbasis+1, nhomo_mon-1);
%Triply_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H10)+abs(H11)+abs(H12)+abs(H20)+abs(H21)+abs(H22))/9
%Triply_degenerate_LUMO_coupling=(abs(L00)+abs(L01)+abs(L02)+abs(L10)+abs(L11)+abs(L12)+abs(L20)+abs(L21)+abs(L22))/9

%H03=H_eff(nhomo_mon+nbasis,   nhomo_mon-3);
%H04=H_eff(nhomo_mon+nbasis,   nhomo_mon-4);
%H13=H_eff(nhomo_mon+nbasis-1, nhomo_mon-3);
%H14=H_eff(nhomo_mon+nbasis-1, nhomo_mon-4);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-3);
%H22=H_eff(nhomo_mon+nbasis-2, nhomo_mon-4);


%Four_fold_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H03)
%				  + abs(H10)+abs(H11)+abs(H12)+abs(H13)
%				  + abs(H20)+abs(H21)+abs(H22)+abs(H23)
%				  + abs(H30)+abs(H31)+abs(H32)+abs(H33))/16
