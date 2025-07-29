function d = pointToSegmentDistance(p, center, orientation, length)
    % Returns minimum distance from point `p` to the line segment defined by center and orientation
    half_vec = 0.5 * length * orientation;
    a = center - half_vec; 
    b = center + half_vec; 

    % Project point onto segment
    ab = b - a;
    ap = p - a;
    t = dot(ap, ab) / dot(ab, ab);
    t = max(0, min(1, t));  % Clamp to segment
    closest = a + t * ab;

    d = norm(p - closest);
end
