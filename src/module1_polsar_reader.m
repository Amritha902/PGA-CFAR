function [SHH, SHV, SVV, inc_angle_map, geo_info] = ...
    module1_polsar_reader(data_path, data_format)
% MODULE1_POLSAR_READER
%   Reads quad-polarimetric SAR data and returns complex scattering
%   matrix channels SHH, SHV, SVV plus incidence angle map.
%
%   Supported formats:
%     'RADARSAT2'  - Requires SHH.bin, SHV.bin, SVV.bin in data_path
%                    and product.xml for metadata
%     'ALOS2'      - CEOS format IMG-HH-*, IMG-HV-*, IMG-VV-*
%     'GAOFEN3'    - GeoTIFF quad-pol (4-band complex)
%     'SYNTHETIC'  - Generates a synthetic scene for code testing

geo_info = struct('lat_ul',0,'lon_ul',0,'lat_lr',0,'lon_lr',0, ...
                  'pixel_spacing_m',10,'format',data_format);

switch upper(data_format)

    % ------------------------------------------------------------------
    case 'RADARSAT2'
        fprintf('  Format: RADARSAT-2 (binary .bin + product.xml)\n');
        % Parse product.xml for scene dimensions
        xml_file = fullfile(data_path, 'product.xml');
        if ~exist(xml_file,'file')
            error('product.xml not found in %s', data_path);
        end
        xml  = xmlread(xml_file);
        nLines = str2double(xml.getElementsByTagName('numberOfLines') ...
                    .item(0).item(0).getData());
        nSamps = str2double(xml.getElementsByTagName('numberOfSamplesPerLine') ...
                    .item(0).item(0).getData());

        SHH = read_rs2_channel(fullfile(data_path,'imagery_HH.bin'), nLines, nSamps);
        SHV = read_rs2_channel(fullfile(data_path,'imagery_HV.bin'), nLines, nSamps);
        SVH = read_rs2_channel(fullfile(data_path,'imagery_VH.bin'), nLines, nSamps);
        SVV = read_rs2_channel(fullfile(data_path,'imagery_VV.bin'), nLines, nSamps);
        % Enforce reciprocity
        SHV = (SHV + SVH) / 2;
        inc_angle_map = load_rs2_incidence(data_path, nLines, nSamps);
        geo_info.format = 'RADARSAT2';

    % ------------------------------------------------------------------
    case 'ALOS2'
        fprintf('  Format: ALOS-2 PALSAR-2 (CEOS)\n');
        d = dir(fullfile(data_path,'IMG-HH-*'));
        if isempty(d), error('No IMG-HH-* file found in %s', data_path); end
        [SHH, nLines, nSamps] = read_ceos_channel(fullfile(data_path, d(1).name));
        d = dir(fullfile(data_path,'IMG-HV-*'));
        SHV = read_ceos_channel(fullfile(data_path, d(1).name), nLines, nSamps);
        d = dir(fullfile(data_path,'IMG-VV-*'));
        SVV = read_ceos_channel(fullfile(data_path, d(1).name), nLines, nSamps);
        inc_angle_map = build_flat_incidence(nLines, nSamps, 35);
        geo_info.format = 'ALOS2';

    % ------------------------------------------------------------------
    case 'GAOFEN3'
        fprintf('  Format: Gaofen-3 (GeoTIFF)\n');
        files = dir(fullfile(data_path,'*.tif'));
        if length(files) < 3
            error('Expected at least 3 GeoTIFF files in %s', data_path);
        end
        SHH = double(imread(fullfile(data_path, files(1).name)));
        SHH = complex(SHH(:,:,1), SHH(:,:,2));
        SHV = double(imread(fullfile(data_path, files(2).name)));
        SHV = complex(SHV(:,:,1), SHV(:,:,2));
        SVV = double(imread(fullfile(data_path, files(3).name)));
        SVV = complex(SVV(:,:,1), SVV(:,:,2));
        [nLines, nSamps] = size(SHH);
        inc_angle_map = build_flat_incidence(nLines, nSamps, 30);
        geo_info.format = 'GAOFEN3';

    % ------------------------------------------------------------------
    case 'SYNTHETIC'
        fprintf('  Format: SYNTHETIC (testing scene)\n');
        [SHH, SHV, SVV, inc_angle_map] = generate_synthetic_scene();
        geo_info.format = 'SYNTHETIC';

    otherwise
        error('Unknown data format: %s', data_format);
end

fprintf('  Data loaded. Size: %d x %d\n', size(SHH,1), size(SHH,2));
end


% =========================================================================
%  RADARSAT-2 helpers
% =========================================================================
function S = read_rs2_channel(fname, nLines, nSamps)
    fid = fopen(fname,'r','b');  % big-endian
    if fid < 0, error('Cannot open: %s', fname); end
    raw = fread(fid, nLines*nSamps*2, 'float32');
    fclose(fid);
    raw = reshape(raw, 2, nSamps, nLines);
    S = squeeze(raw(1,:,:) + 1j*raw(2,:,:))';
end

function inc = load_rs2_incidence(data_path, nLines, nSamps)
    % Try to find lutSigma.xml or incidenceAngle
    xml_file = fullfile(data_path, 'product.xml');
    try
        xml = xmlread(xml_file);
        nodes = xml.getElementsByTagName('incidenceAngle');
        if nodes.getLength() > 0
            vals = str2num(char(nodes.item(0).item(0).getData())); %#ok
            inc  = repmat(linspace(vals(1),vals(end),nSamps), nLines, 1);
            return;
        end
    catch
    end
    inc = build_flat_incidence(nLines, nSamps, 35);
end


% =========================================================================
%  ALOS-2 CEOS helper (simplified — assumes SLC format)
% =========================================================================
function [S, nLines, nSamps] = read_ceos_channel(fname, nL, nS)
    fid = fopen(fname,'r','b');
    if fid < 0, error('Cannot open: %s', fname); end
    % Skip 720-byte file descriptor
    fseek(fid, 720, 'bof');
    % Each record: 412-byte prefix + nSamps*2*2 bytes (int16 I/Q)
    if nargin < 2
        % Read first record prefix to infer dimensions
        fseek(fid, 720+180, 'bof');
        nSamps = fread(fid,1,'uint32');
        fseek(fid, 720, 'bof');
        raw_all = fread(fid, Inf, 'int16');
        total = length(raw_all);
        rec_bytes = (412 + nSamps*4);
        nLines = floor(total*2 / rec_bytes);
    else
        nLines = nL; nSamps = nS;
    end
    fseek(fid, 720, 'bof');
    S = complex(zeros(nLines, nSamps,'single'));
    for l = 1:nLines
        fseek(fid, 412, 'cof');  % skip record header
        raw = fread(fid, nSamps*2, 'int16');
        if length(raw) < nSamps*2, break; end
        S(l,:) = single(raw(1:2:end)) + 1j*single(raw(2:2:end));
    end
    fclose(fid);
end


% =========================================================================
%  Utility: flat incidence angle ramp
% =========================================================================
function inc = build_flat_incidence(nRows, nCols, center_deg)
    inc = repmat(linspace(center_deg-5, center_deg+5, nCols), nRows, 1);
end


% =========================================================================
%  SYNTHETIC SCENE GENERATOR
%  Creates a 512x512 test scene with:
%   - Calm ocean (low PSSI)
%   - Rough ocean patch (high PSSI)
%   - 8 ship-like targets
% =========================================================================
function [SHH, SHV, SVV, inc_angle_map] = generate_synthetic_scene()

    fprintf('    Generating synthetic 512x512 quad-pol scene...\n');
    N = 512;
    rng(42);

    % Incidence angle ramp (25–45 degrees across range)
    inc_angle_map = repmat(linspace(25,45,N), N, 1);
    theta_r = deg2rad(inc_angle_map);

    % ---- Sea clutter parameters ----
    % Bragg resonance model: sigma_HH / sigma_VV varies with incidence
    eps_r  = 70 - 1j*40;   % Seawater dielectric constant at C-band
    cos_t  = cos(theta_r);
    sin_t  = sin(theta_r);
    eps_cos = sqrt(eps_r - sin_t.^2);

    % Fresnel coefficients
    R_HH = (cos_t - eps_cos) ./ (cos_t + eps_cos);
    R_VV = (eps_r.*cos_t - eps_cos) ./ (eps_r.*cos_t + eps_cos);
    Bpred = abs(R_HH ./ R_VV).^2;   % Predicted HH/VV ratio

    % ---- Calm ocean region (upper 70%) ----
    calm_mask = true(N,N);
    calm_mask(round(0.7*N)+1:end, :) = false;

    % Amplitude from K-distribution (calm sea)
    k_shape_calm = 5.0;
    sigma0_HH_calm = 0.01;
    amp_HH_calm = sqrt(sigma0_HH_calm) * krnd(N, N, k_shape_calm);
    amp_VV_calm = amp_HH_calm .* sqrt(Bpred);
    amp_HV_calm = amp_HH_calm * 0.15;  % Cross-pol 15 dB below co-pol

    % Rough ocean patch (lower 30%)
    k_shape_rough = 1.5;
    sigma0_HH_rough = 0.05;
    amp_HH_rough = sqrt(sigma0_HH_rough) * krnd(N, N, k_shape_rough);
    amp_VV_rough = amp_HH_rough .* sqrt(Bpred) * 0.7;  % Ratio deviates
    amp_HV_rough = amp_HH_rough * 0.30;

    % Blend calm/rough
    amp_HH = amp_HH_calm; amp_HH(~calm_mask) = amp_HH_rough(~calm_mask);
    amp_VV = amp_VV_calm; amp_VV(~calm_mask) = amp_VV_rough(~calm_mask);
    amp_HV = amp_HV_calm; amp_HV(~calm_mask) = amp_HV_rough(~calm_mask);

    % Random phases — calm: CPD ~ 0; rough: CPD deviates
    phi_HH = 2*pi*rand(N,N);
    phi_VV_calm  = phi_HH + (pi/18)*randn(N,N);   % ~0 deg CPD
    phi_VV_rough = phi_HH + pi + (pi/6)*randn(N,N); % ~180 deg CPD region
    phi_VV = phi_VV_calm; phi_VV(~calm_mask) = phi_VV_rough(~calm_mask);
    phi_HV = 2*pi*rand(N,N);

    SHH = amp_HH .* exp(1j*phi_HH);
    SVV = amp_VV .* exp(1j*phi_VV);
    SHV = amp_HV .* exp(1j*phi_HV);

    % ---- Insert 8 ships ----
    % Ship positions (row, col), approximate size in pixels
    ships = [100,  80, 6, 4;   % [row, col, len, wid]
             130, 200, 8, 3;
             200, 350, 5, 4;
              80, 430, 7, 3;
             300, 100, 6, 5;
             280, 300, 9, 4;
             350, 450, 5, 3;
             420, 200, 7, 4];

    for s = 1:size(ships,1)
        r  = ships(s,1); c  = ships(s,2);
        rl = ships(s,3); rw = ships(s,4);
        rr = max(1,r-rl):min(N,r+rl);
        cc = max(1,c-rw):min(N,c+rw);
        % Double-bounce signature: high |b|^2, low H, CPD ~ 180
        amp_ship   = 0.8 + 0.3*rand();
        phase_ship = pi + (pi/12)*randn(length(rr), length(cc));
        SHH(rr,cc) = amp_ship * exp(1j*phi_HH(rr,cc));
        SVV(rr,cc) = amp_ship * exp(1j*(phi_HH(rr,cc) + phase_ship));
        SHV(rr,cc) = 0.02 * amp_ship * exp(1j*phi_HV(rr,cc));
    end

    fprintf('    Synthetic scene ready. Ships planted at %d locations.\n', ...
        size(ships,1));
end


function x = krnd(M, N, nu)
% K-distributed random amplitudes (texture * speckle)
    texture = gamrnd(nu, 1/nu, M, N);
    speckle = abs(randn(M,N) + 1j*randn(M,N)) / sqrt(2);
    x = sqrt(texture) .* speckle;
end
