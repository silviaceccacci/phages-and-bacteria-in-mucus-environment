function isAttached = check_PhageToRodAttachment(phage_pos, bact_center, lB, d_enc1)
    
        %Compute endpoints of the bacterium rod
        half_L = 0.5 * lB;
        a = bact_center - half_L;
        b = bact_center + half_L;

        %Compute vector from a to b and from a to phage
        v = b - a;
        w = phage_pos - a;

        %Project w onto v and clamp to [0,1]
        t = dot(w,v) / dot(v,v);
        t = max(0, min(1, t));

        %Closest point on the rod of the phage
        closest_point = a + t*v;

        %Compute distance from phage center to closest point
        d = norm(phage_pos - closest_point);
       
        isAttached = (d < d_enc1); %Return true if within threshold

end