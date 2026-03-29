function [det_map, thr_map] = run_h_cfar(Span_map, H_map, params)
%RUN_H_CFAR Entropy-adaptive k, fixed square window.
k_map = params.kmin + (params.kmax - params.kmin) .* double(H_map);
r0    = params.r0;
guard = max(1, floor(r0/3));

[nR,nC] = size(Span_map);
thr_map = zeros(nR,nC,'single');
det_map = false(nR,nC);
pad     = r0+1;
Sp      = padarray(double(Span_map),[pad pad],'replicate','both');

hw = r0; gw = guard;
outer_mask = true(2*hw+1,2*hw+1);
outer_mask(hw+1-gw:hw+1+gw, hw+1-gw:hw+1+gw) = false;
outer_mask(hw+1,hw+1) = false;
bg_idx = find(outer_mask);
[bg_dx,bg_dy] = ind2sub(size(outer_mask),bg_idx);
bg_dx = bg_dx-(hw+1); bg_dy = bg_dy-(hw+1);

for ri = 1:nR
    for ci = 1:nC
        rp=ri+pad; cp=ci+pad;
        bgr=rp+bg_dx; bgc=cp+bg_dy;
        bg_vals=Sp(sub2ind(size(Sp),bgr,bgc));
        mu=mean(bg_vals); sig=std(bg_vals);
        thr=mu+k_map(ri,ci)*sig;
        thr_map(ri,ci)=thr;
        det_map(ri,ci)=Span_map(ri,ci)>thr;
    end
end
end
