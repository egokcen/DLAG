function B = pairwise_coeffs_pcca(in, out, params)
%
% B = pairwise_coeffs_pcca(in, out, params)
%
% Description: Compute the reduced-rank regression coefficient matrix B,
%              between the group indexed by 'out' and the group indexed
%              by 'in', given pCCA model params.
%
% Arguments:
%     in     -- int; Index of group used for prediction.
%     out    -- int; Index of group to be predicted.
%     params -- learned PCCA parameteres (structure with fields Cs, Rs, ds)
%
% Outputs:
%     B -- (out_dim x in_dim) array; pairwise reduced-rank regression
%          coefficient matrix.
%
% Author: Evren Gokcen

    Cin = params.Cs{in};
    Rin = params.Rs{in};
    din = params.ds{in};
    
    Cout = params.Cs{out};
    dout = params.ds{out};
    
    B = Cout * Cin' / (Cin * Cin' + Rin);

end