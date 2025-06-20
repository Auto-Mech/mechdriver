*gen40

  title
    ${title}
  end

  curvature compute
  gradder noextra
  freqsource hessian
  hessform full
  anharm none
  maxlpts 3
  npts ${npts}
  hessunit bohr
  geomunit bohr
  gradunit bohr
  print debugging


${dat_str}
