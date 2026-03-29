function ais_data = load_ais_csv(ais_file, geo_info, nRows, nCols)
%LOAD_AIS_CSV Reads AIS CSV and converts lat/lon to pixel coordinates.
%  CSV format: columns Lat, Lon, MMSI (with header row)

T = readtable(ais_file);

ais_data.lat  = T.Lat;
ais_data.lon  = T.Lon;
ais_data.mmsi = T.MMSI;

lat_range = geo_info.lat_ul - geo_info.lat_lr;
lon_range = geo_info.lon_lr - geo_info.lon_ul;

ais_data.row_px = round((geo_info.lat_ul - T.Lat) ./ (lat_range + eps) * nRows);
ais_data.col_px = round((T.Lon - geo_info.lon_ul) ./ (lon_range + eps) * nCols);

valid = ais_data.row_px >= 1 & ais_data.row_px <= nRows & ...
        ais_data.col_px >= 1 & ais_data.col_px <= nCols;

ais_data.row_px = ais_data.row_px(valid);
ais_data.col_px = ais_data.col_px(valid);
ais_data.lat    = ais_data.lat(valid);
ais_data.lon    = ais_data.lon(valid);

fprintf('  AIS: %d total ships, %d within scene bounds.\n', height(T), sum(valid));
end
