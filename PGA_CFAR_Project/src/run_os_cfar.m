function [det_map, thr_map] = run_os_cfar(Span_map, params)
%RUN_OS_CFAR Ordered-Statistics CFAR (75th percentile).
k_fixed  = (params.kmin + params.kmax) / 2;
r0       = params.r0;
guard    = max(1, floor(r0/3));
pct_rank = 0.75;

[nR, nC] = size(Span_map);
thr_map  = zeros(nR, nC, 'single');
det_map  = false(nR, nC);
pad      = r0 + 1;
Sp       = padarray(double(Span_map), [pad pad], 'replicate', 'both');

hw = r0; gw = guard;
[~,~] = meshgrid(-hw:hw,-hw:hw);
outer_mask = true(2*hw+1,2*hw+1);
outer_mask(-gw+hw+1:gw+hw+1,-gw+hw+1:gw+hw+1)=false;
outer_mask(hw+1,hw+1)=false;
bg_idx = find(outer_mask);
[bg_dx,bg_dy] = ind2sub(size(outer_mask),bg_idx);
bg_dx = bg_dx-(hw+1); bg_dy = bg_dy-(hw+1);

for ri = 1:nR
    for ci = 1:nC
        rp = ri+pad; cp = ci+pad;
        bgr = rp+bg_dx; bgc = cp+bg_dy;
        bg_vals = sort(Sp(sub2ind(size(Sp),bgr,bgc)));
        n   = length(bg_vals);
        kth = bg_vals(min(n, max(1,round(pct_rank*n))));
        thr = k_fixed * kth;
        thr_map(ri,ci) = thr;
        det_map(ri,ci) = Span_map(ri,ci) > thr;
    end
end
end


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
outer_mask(-gw+hw+1:gw+hw+1,-gw+hw+1:gw+hw+1)=false;
outer_mask(hw+1,hw+1)=false;
bg_idx = find(outer_mask);
[bg_dx,bg_dy] = ind2sub(size(outer_mask),bg_idx);
bg_dx=bg_dx-(hw+1); bg_dy=bg_dy-(hw+1);

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


function [det_map, thr_map] = run_pga_cfar_threshold_only(Span_map, PSSI_map, ...
                                    alpha_map, H_map, params)
%RUN_PGA_CFAR_THRESHOLD_ONLY PSSI-adaptive k + pre-screen, fixed square window.
k_map     = params.kmin + (params.kmax-params.kmin) .* double(PSSI_map);
prescreen = (alpha_map > params.alpha_thresh) & (H_map < params.H_thresh);
r0        = params.r0;
guard     = max(1,floor(r0/3));

[nR,nC]  = size(Span_map);
thr_map  = zeros(nR,nC,'single');
det_map  = false(nR,nC);
pad      = r0+1;
Sp       = padarray(double(Span_map),[pad pad],'replicate','both');

hw=r0; gw=guard;
outer_mask=true(2*hw+1,2*hw+1);
outer_mask(-gw+hw+1:gw+hw+1,-gw+hw+1:gw+hw+1)=false;
outer_mask(hw+1,hw+1)=false;
bg_idx=find(outer_mask);
[bg_dx,bg_dy]=ind2sub(size(outer_mask),bg_idx);
bg_dx=bg_dx-(hw+1); bg_dy=bg_dy-(hw+1);

for ri=1:nR
    for ci=1:nC
        if ~prescreen(ri,ci), continue; end
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
