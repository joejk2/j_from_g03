disp ("Assumes the presence of the following files: MOs_1.txt, MOs_2.txt, MOs_pair.txt, S_1.txt, S_2.txt, S_pair.txt, Evls_pair.txt")
disp ("!!! Octave counts from 1 !!!")

%nbasis   =input("Number of MOs (of monomer)");		%... of monomer
%nhomo_mon=input("HOMO (of monomer)");		%... of monomer (!! Octave counts from 1 !!)
%ethylene
nbasis=46
nhomo_mon=8

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
should_be_zero=H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis-1)
should_be_zero=H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis-2)
should_be_zero=H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis-3)
should_be_zero=H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis-4)
not_zero=H_eff(nhomo_mon-1,nhomo_mon-1)
not_zero=H_eff(nhomo_mon-2,nhomo_mon-2)
not_zero=H_eff(nhomo_mon-3,nhomo_mon-3)
not_zero=H_eff(nhomo_mon-4,nhomo_mon-4)
H00=H_eff(nhomo_mon+nbasis,   nhomo_mon  );
E1_E2=(H_eff(nhomo_mon,nhomo_mon)-H_eff(nhomo_mon+nbasis,nhomo_mon+nbasis))
half_HOMO_diff=(Evls(2*nhomo_mon,2*nhomo_mon)-Evls(2*nhomo_mon-1,2*nhomo_mon-1))/2
Non_degenerate_HOMO_coupling=H_eff(nhomo_mon+nbasis,   nhomo_mon  )
Non_degenerate_LUMO_coupling=H_eff(nhomo_mon+1+nbasis, nhomo_mon+1)
%L00=H_eff(nhomo_mon+nbasis+1, nhomo_mon+1);
%Non_degenerate_LUMO_coupling=abs(L00)

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
%Five_fold_degenerate_HOMO_coupling=(abs(H00)+abs(H01)+abs(H02)+abs(H03)+abs(H04)
%				  + abs(H10)+abs(H11)+abs(H12)+abs(H13)+abs(H14)
%				  + abs(H20)+abs(H21)+abs(H22)+abs(H23)+abs(H24)
%				  + abs(H30)+abs(H31)+abs(H32)+abs(H33)+abs(H34)
%				  + abs(H40)+abs(H41)+abs(H42)+abs(H43)+abs(H44))/25
