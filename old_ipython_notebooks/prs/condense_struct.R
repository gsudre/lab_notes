regional = read.csv('~/data/prs/REGIONAL_ANALYSES_FREESURFER.csv')

condense_sublobar = function(ldata) {
  sublobar_regions = unique(regional$subloar)
  # for each region, add or average the columns corresponding to the ROIs
  area_data = c()
  hdr = c()
  volume_data = c()
  thickness_data = c()
  aparc_data = c()
  aparc_hdr = c()
  for (s in sublobar_regions) {
    if (s != "") {
      idx = which(regional$subloar == s)
      area = 0
      volume = 0
      thickness = 0
      aparc = 0
      for (i in idx) {
        cnt = 0
        roi_name = sub('_thickness', '', regional[i, 1])
        if (length(grep('rh_|lh_', roi_name)) > 0) {
          var_name = quote(sprintf('%s_area', roi_name))
          area = area + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_volume', roi_name))
          volume = volume + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_thickness', roi_name))
          thickness = thickness + ldata[, eval(var_name)]
          cnt = cnt + 1
        }
        else {
          var_name = quote(roi_name)
          aparc = aparc + ldata[, eval(var_name)]
        }
      }
      if (sum(area) > 0) {
        hdr = c(hdr, s)
        area_data = cbind(area_data, area)
        volume_data = cbind(volume_data, volume)
        thickness_data = cbind(thickness_data, thickness / cnt)
      }
      else {
        aparc_hdr = c(aparc_hdr, s)
        aparc_data = cbind(aparc_data, aparc)
      }
    }
  }
  colnames(area_data) = lapply(hdr, function(d) sprintf('%s_area', d))
  colnames(volume_data) = lapply(hdr, function(d) sprintf('%s_volume', d))
  colnames(thickness_data) = lapply(hdr, function(d) sprintf('%s_thickness', d))
  colnames(aparc_data) = aparc_hdr
  return(cbind(area_data, volume_data, thickness_data, aparc_data))
}

condense_lobar = function(ldata) {
  lobar_regions = unique(regional$lobar)
  # for each region, add or average the columns corresponding to the ROIs
  areaLobar_data = c()
  hdr = c()
  volumeLobar_data = c()
  thicknessLobar_data = c()
  aparcLobar_data = c()
  aparc_hdr = c()
  for (s in lobar_regions) {
    if (s != "") {
      idx = which(regional$lobar == s)
      area = 0
      volume = 0
      thickness = 0
      aparc = 0
      for (i in idx) {
        cnt = 0
        roi_name = sub('_thickness', '', regional[i, 1])
        if (length(grep('rh_|lh_', roi_name)) > 0) {
          var_name = quote(sprintf('%s_area', roi_name))
          area = area + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_volume', roi_name))
          volume = volume + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_thickness', roi_name))
          thickness = thickness + ldata[, eval(var_name)]
          cnt = cnt + 1
        }
        else {
          var_name = quote(roi_name)
          aparc = aparc + ldata[, eval(var_name)]
        }
      }
      if (sum(area) > 0) {
        hdr = c(hdr, s)
        areaLobar_data = cbind(areaLobar_data, area)
        volumeLobar_data = cbind(volumeLobar_data, volume)
        thicknessLobar_data = cbind(thicknessLobar_data, thickness / cnt)
      }
      else {
        aparc_hdr = c(aparc_hdr, s)
        aparcLobar_data = cbind(aparcLobar_data, aparc)
      }
    }
  }
  colnames(areaLobar_data) = lapply(hdr, function(d) sprintf('%s_area', d))
  colnames(volumeLobar_data) = lapply(hdr, function(d) sprintf('%s_volume', d))
  colnames(thicknessLobar_data) = lapply(hdr, function(d) sprintf('%s_thickness', d))
  colnames(aparcLobar_data) = aparc_hdr
  return(cbind(areaLobar_data, volumeLobar_data,
               thicknessLobar_data, aparcLobar_data))
}

condense_theory = function(ldata) {
  theory_regions = unique(regional$ROI_theory_driven)
  # for each region, add or average the columns corresponding to the ROIs
  areaTheory_data = c()
  hdr = c()
  volumeTheory_data = c()
  thicknessTheory_data = c()
  aparcTheory_data = c()
  aparc_hdr = c()
  for (s in theory_regions) {
    if (s != "") {
      idx = which(regional$ROI_theory_driven == s)
      area = 0
      volume = 0
      thickness = 0
      aparc = 0
      for (i in idx) {
        cnt = 0
        roi_name = sub('_thickness', '', regional[i, 1])
        if (length(grep('rh_|lh_', roi_name)) > 0) {
          var_name = quote(sprintf('%s_area', roi_name))
          area = area + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_volume', roi_name))
          volume = volume + ldata[, eval(var_name)]
          var_name = quote(sprintf('%s_thickness', roi_name))
          thickness = thickness + ldata[, eval(var_name)]
          cnt = cnt + 1
        }
        else {
          var_name = quote(roi_name)
          aparc = aparc + ldata[, eval(var_name)]
        }
      }
      if (sum(area) > 0) {
        hdr = c(hdr, s)
        areaTheory_data = cbind(areaTheory_data, area)
        volumeTheory_data = cbind(volumeTheory_data, volume)
        thicknessTheory_data = cbind(thicknessTheory_data, thickness / cnt)
      }
      else {
        aparc_hdr = c(aparc_hdr, s)
        aparcTheory_data = cbind(aparcTheory_data, aparc)
      }
    }
  }
  colnames(areaTheory_data) = lapply(hdr, function(d) sprintf('%s_area', d))
  colnames(volumeTheory_data) = lapply(hdr, function(d) sprintf('%s_volume', d))
  colnames(thicknessTheory_data) = lapply(hdr, function(d) sprintf('%s_thickness', d))
  colnames(aparcTheory_data) = aparc_hdr
  return(cbind(areaTheory_data, volumeTheory_data,
               thicknessTheory_data, aparcTheory_data))
}
