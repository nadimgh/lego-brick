function [ u, x, LLR ] = polarDecodeSC( L, info_set )
% Polar successive cancellation decoder
%   L = 1xN vector of log-likelihood ratios for the received sequence 
%   info_set = 1xN vector with ones in the positions of non-frozen indices and zeros otherwise
%   u = 1xN vector for the input sequence to the polar encoder
%   x = 1xN vector for the input sequence to the channel (output of the polar encoder)
%   LLR = 1xN vector for the log-likelihood ratio of polarized channels

N = length(L);
if (N==1)
    LLR = L;
    if info_set == 0
        % If frozen bit, set output to 0
        x = 0; u = x;
    else
        % If non-frozen bit, check log-likelihood
        if (L<0)
            x = 1; u = x;
        else
            x = 0; u = x;
        end
    end
else
    info_set_1 = info_set(1:N/2);
    info_set_2 = info_set(N/2+1:end);
    L1 = f(L(1:N/2), L(N/2+1:end));
    [uhat1, vhat1, LLR1] = polarDecodeSC(L1, info_set_1);
    L2 = g(L(1:N/2), L(N/2+1:end), vhat1);
    [uhat2, vhat2, LLR2] = polarDecodeSC(L2, info_set_2);
    LLR = [LLR1 LLR2];
    u = [uhat1 uhat2];
    x = zeros(1,N);
    x(1:N/2) = mod(vhat1+vhat2, 2);
    x(N/2+1:end) = vhat2;
end
return

function z = f(x,y)
z = 2*atanh(tanh(x/2).*tanh(y/2));
return

% function z = f2(x,y)
% xabs = abs(x);
% yabs = abs(y);
% maxi = max(xabs,yabs);
% mini = min(xabs,yabs);
% z = sign(x).*sign(y).*(mini+log1p(exp(-(xabs+yabs)))-log1p(exp(-(maxi-mini))));
% return

function z = g(x,y,u1)
z = (1-2*u1).*x + y;
return
