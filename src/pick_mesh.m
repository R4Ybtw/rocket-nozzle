function [mesh, P_exit, P_wall, P_axis] = pick_mesh(res, cas, preferred)
%PICK_MESH  Profiles of one case from the preferred mesh, else the finest available.
%   [mesh, P_exit, P_wall, P_axis] = pick_mesh(res, cas, preferred)
    order = {preferred, 'medium', 'fine', 'coarse'};
    for k = 1:numel(order)
        mesh = order{k};
        if isempty(mesh), continue; end
        P_exit = read_profile(fullfile(res, sprintf('%s_%s_exit.csv', cas, mesh)));
        P_wall = read_profile(fullfile(res, sprintf('%s_%s_wall.csv', cas, mesh)));
        P_axis = read_profile(fullfile(res, sprintf('%s_%s_axis.csv', cas, mesh)));
        if P_exit.n + P_wall.n + P_axis.n > 0, return; end
    end
    mesh = '–';
end
