clear

%% Test without rotors
disp('Regular revolutes')

model_reg = Arm6LinkModel();

q = rand(model_reg.NQ,1);
q = model_reg.normalizeConfVec(q);
qd = rand(model_reg.NV,1);
qdd = rand(model_reg.NV,1);

[tau_reg,out_reg] = ID(model_reg,q,qd,qdd);
[tau2_reg,out2_reg] = ID_Lee(model_reg,q,qd,qdd);
etau = norm(tau_reg-tau2_reg)

Y_reg = RegressorClassical(model_reg, q,qd,qdd);
% Y2_reg = RegressorWithDerivs_Lee(model_reg, q,qd,qdd);

% eY = norm( tau - Y*a)
% eY = norm(Y_reg - Y2_reg)


%% Test with rotors
disp('Revolutes with rotors')

model_rot = Arm6LinkRotorModel();

q = rand(model_rot.NQ,1);
q = model_rot.normalizeConfVec(q);
qd = rand(model_rot.NV,1);
qdd = rand(model_rot.NV,1);

[tau_rot,out_rot] = ID(model_rot,q,qd,qdd);
[tau2_rot,out2_rot] = ID_Lee(model_rot,q,qd,qdd);
etau = norm(tau_rot-tau2_rot)

Y_rot = RegressorClassical(model_rot, q,qd,qdd);
% Y2_rot = RegressorWithDerivs_Lee(model_rot, q,qd,qdd);

% eY = norm( tau - Y*a)
% eY = norm(Y_rot - Y2_rot)

%% Test with absolute triplet
disp('Absolute triplet')
model_abs = Arm6LinkAbsModel();

q = rand(model_abs.NQ,1);
q = model_abs.normalizeConfVec(q);
qd = rand(model_abs.NV,1);
qdd = rand(model_abs.NV,1);

[tau_abs,out_abs] = ID(model_abs,q,qd,qdd);
[tau2_abs,out2_abs] = ID_Lee(model_abs,q,qd,qdd);
etau = norm(tau_abs-tau2_abs)

Y_abs = RegressorClassical(model_abs, q,qd,qdd);
% Y2_abs = RegressorWithDerivs_Lee(model_abs, q,qd,qdd);

% eY = norm( tau - Y*a)
% eY = norm(Y_abs - Y2_abs)

%% Testing ID derivatives

% Fourier trajectory parameters
num_sample       = 100;     % number of samples of the trajectory
horizon          = 10;      % trajectory horizon
base_frequency   = pi*0.3;  % Fourier trajectory base frequency
sample_times      = linspace(0,horizon,num_sample);
sample_interval  = horizon/num_sample;              % dt for numerical integration
m         = 5;                       % num of trajectory coefs per joints
p_initial = rand(2*m,6) - 0.5;  % initial p
% params is ordered [a_1, b_1, a_2, b_2, ..., a_m, b_m] for Fourier Series: a_k*sin(kwt) + b_k*cos(kwt)

% Fourier trajectory generation with parameter p_initial
[q, qd, qdd]    = makeFourier(p_initial, base_frequency, sample_times);

[dq, dqd, dqdd] = getFourierDerivative(p_initial, base_frequency, sample_times);
% for each joint i and each time j, taking the derivative w.r.t all 2m
% parameters...where does extra i index come from?

test_index = 23;

qtest = q(:,test_index);
qdtest = qd(:,test_index);
qddtest = qdd(:,test_index);

dqtest = dq(:,:,test_index);
dqdtest = dqd(:,:,test_index);
dqddtest = dqdd(:,:,test_index);

model_to_use = model_abs;

a = [];
a(1:10,1) = inertiaMatToVec( model_to_use.I_RB{1} );
for i = 2:model_to_use.N_RB
    a(end+1 : end+10) = inertiaMatToVec(model_to_use.I_RB{i});
end

[tau,out] = ID(model_to_use,qtest,qdtest,qddtest);
[tau2,out] = ID_Lee(model_to_use,qtest,qdtest,qddtest);

[dtau, dV, dVd] = ID_derivatives_Lee( model_to_use, qtest, qdtest, qddtest, ...
                        dqtest, dqdtest, dqddtest, out.f );

[Y_out, W_out] = RegressorClassical(model_to_use, qtest,qdtest,qddtest);
[Y2_out, dY2_out, W2_out, dW2_out] = RegressorWithDerivs_Lee(model_to_use, qtest,qdtest,qddtest, ...
                                    dqtest,dqdtest,dqddtest);

eW = norm(W_out - W2_out)
eY = norm(Y_out - Y2_out)

eY1 = norm( tau - Y_out*a)
eY2 = norm( tau2 - Y2_out*a)