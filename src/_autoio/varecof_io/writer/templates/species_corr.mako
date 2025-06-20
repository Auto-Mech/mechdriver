
      real*8 function ${species_name}_corr(natoms,x,rparm,iparm)
c     natoms = number of atoms (not used)
c     x = cartesian coordinates
c     rparm (not used)
c     iparm(1) = specifies which potential correction to use
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(${npot})
      dimension rinp(${npot_terms})
      dimension ${dv_defs}
      dimension dv20(${npot_terms})
c     rinp = A-B distance, where A,B are the bond forming atoms
${rvals}
      data nrin / ${npot_terms} /
c     dvN = energy value for correction potential N
% if pot_labels != '':
${pot_labels}
% endif
      data dvp1,dvpn / 1.0d40,1.0d40 /
${dv_vals}
## Append further info that is needed for atom labeling
c     min and max values along the rAB distance grid
      data rABmin,rABmax / ${rmin},${rmax} /

c     Set the distance of the two atoms forming the bond in the reaction
${bond_distance_string}

## Set strings that modify the potential for certain distances
% if restrict_distance_strings != '':
c     Set the distances to other atoms we must check to alter the potential
${restrict_distance_strings}
% endif
c     Set mult factor for potential to exponential form when sampling dist outside bounds
${delmlt_string}

## Append the lines for declaring the spline functions
c     Perform spline fits for each of the correction potential grids
      ipot = iparm(1)
${spline_strings}

c     Set the value of the correction potential
      ${species_name}_corr = ${species_name}_corr*delmlt/627.5095

      return
      end

