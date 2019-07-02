% 論文「動加速度環境下における 姿勢推定アルゴリズムの研究」
% EKFアルゴリズムの数式計算

clear

syms q1 q2 q3 q4 ob1 ob2 ob3 om1 om2 om3 b1 b2 b3 D a1 a2 a3 m1 m2 m3
q = [q1; q2; q3; q4];%姿勢クォータニオン
ob = [ob1; ob2; ob3];% 角速度バイアス
x = [q; ob];%状態ベクトル
om = [om1; om2; om3];% 角速度
a = [a1; a2; a3];%初期姿勢での加速度
m = [m1; m2; m3];%初期姿勢での磁気


qmat = @(p) [p(1) -p(2) -p(3) -p(4); p(2) p(1) -p(4) p(3); p(3) p(4) p(1) -p(2); p(4) -p(3) p(2) p(1)] % クォータニオン積のための(9)式の行列をつくる 
qmulti = @(q, p) qmat(q)*p %(9)クォータニオンの積
qadj = @(p) [p(1); -p(2); -p(3); -p(4)];%共役クォータニオンをつくる 式(10)
trans = @(r, q) qmulti(qmulti(qadj(q),[0; r]), q)%論文(11)式：r-frameの代数ベクトルrをb-frame上代数ベクトルに座標変換（b-frameの姿勢：クォータニオンq）

% f = [1/2 * qmulti(q, [0; wm]) - 1/2 * qmulti(q, [0; ob]);...
%  diag([-b1 -b2 -b3]) * ob]% (27)式

% 連続時間ダイナミクス　dot{x} = f(x) + G * w
f = [1/2 * qmulti(q, [0; om - ob]); diag([-b1 -b2 -b3]) * ob]% (27)式
G = [zeros(4,3);eye(3)]

% 離散時間ダイナミクス　x_{t+1} = x_{t} + f(x_{t}) * Delta + G * Delta  * w
ft = x + f * D % (30)(31)式
Gt = G*D

F = jacobian(ft, x)
 


h = [qmulti(qadj(q),qmulti([0; a], q)); qmulti(qadj(q),qmulti([0; m], q))]; %(35)式

H = jacobian(h, x)


replace_C(ft)
replace_C(h)
replace_C(F)
replace_C(H)

function ft = replace_C(ft)
% Cプログラム用に文字置換

%置換リスト
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



