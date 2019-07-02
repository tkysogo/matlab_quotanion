% �_���u�������x�����ɂ����� �p������A���S���Y���̌����v
% EKF�A���S���Y���̐����v�Z

clear

syms q1 q2 q3 q4 ob1 ob2 ob3 om1 om2 om3 b1 b2 b3 D a1 a2 a3 m1 m2 m3
q = [q1; q2; q3; q4];%�p���N�H�[�^�j�I��
ob = [ob1; ob2; ob3];% �p���x�o�C�A�X
x = [q; ob];%��ԃx�N�g��
om = [om1; om2; om3];% �p���x
a = [a1; a2; a3];%�����p���ł̉����x
m = [m1; m2; m3];%�����p���ł̎��C


qmat = @(p) [p(1) -p(2) -p(3) -p(4); p(2) p(1) -p(4) p(3); p(3) p(4) p(1) -p(2); p(4) -p(3) p(2) p(1)] % �N�H�[�^�j�I���ς̂��߂�(9)���̍s������� 
qmulti = @(q, p) qmat(q)*p %(9)�N�H�[�^�j�I���̐�
qadj = @(p) [p(1); -p(2); -p(3); -p(4)];%�����N�H�[�^�j�I�������� ��(10)
trans = @(r, q) qmulti(qmulti(qadj(q),[0; r]), q)%�_��(11)���Fr-frame�̑㐔�x�N�g��r��b-frame��㐔�x�N�g���ɍ��W�ϊ��ib-frame�̎p���F�N�H�[�^�j�I��q�j

% f = [1/2 * qmulti(q, [0; wm]) - 1/2 * qmulti(q, [0; ob]);...
%  diag([-b1 -b2 -b3]) * ob]% (27)��

% �A�����ԃ_�C�i�~�N�X�@dot{x} = f(x) + G * w
f = [1/2 * qmulti(q, [0; om - ob]); diag([-b1 -b2 -b3]) * ob]% (27)��
G = [zeros(4,3);eye(3)]

% ���U���ԃ_�C�i�~�N�X�@x_{t+1} = x_{t} + f(x_{t}) * Delta + G * Delta  * w
ft = x + f * D % (30)(31)��
Gt = G*D

F = jacobian(ft, x)
 


h = [qmulti(qadj(q),qmulti([0; a], q)); qmulti(qadj(q),qmulti([0; m], q))]; %(35)��

H = jacobian(h, x)


replace_C(ft)
replace_C(h)
replace_C(F)
replace_C(H)

function ft = replace_C(ft)
% C�v���O�����p�ɕ����u��

%�u�����X�g
rep_lst =[...
    "q1",   "q2",   "q3",   "q4",   "ob1",  "ob2",     "ob3", "om1",    "om2",     "om3", "b1",     "b2", "b3",     "D", "a1", "a2", "a3", "m1", "m2", "m3";
    "q[0]", "q[1]", "q[2]", "q[3]", "ob[0]", "ob[1]", "ob[2]", "om[0]", "om[1]", "om[2]", "b[0]", "b[1]", "b[2]", "D", "a[0]", "a[1]", "a[2]", "m[0]", "m[1]", "m[2]"];

[row column]=size(ft);
ft=string(ft);
for k = 1:row
    for l=1:column
        ft(k,l) = string(ft(k,l));
        for j=1:length(rep_lst)
            ft(k,l)=strrep(ft(k,l), rep_lst(1,j), rep_lst(2,j));
        end
    end
end
end

% fts=sym2cell(ft)
% fts{1,1}=strrep(fts{1,1}, "q1", "q[1]");
% fts{1,1}=strrep(fts{1,1}, "q2", "q[2]");
% 
% fts{1,1}



