 ;***********************************************************************************************************************
; Pre-procesador químico del OLE2-CMC, generador del archivo de emisiones necesario para ejecutar el modelo WRF-CHEM.
; Copyright: Muñoz-Lomas-Vara
;***********************************************************************************************************************
; Instrucciones de uso: Se recomienda modificar apenas dos variables, 'npts' que representa el numero total de puntos 
; emisores y 'data' que controla la existencia de los archivos de entrada relacionados a cinco especies quimicas (NO, 
; CO, SO2, VOC y PM), siendo los VOC y PM posteriormente fraccionados (experimentalmente o de acuerdo a la literatura, 
; ya que por lo general se suelen medir hidrocarbonetos y material particulado de forma compacta) en 16 y 11 
; subespecies, completando asi 30 especies, si se quiere ejecutar usando el mecanismo quimico RADM2 por ejemplo. Usar 
; "si" y por ejemplo "no" segun sea el caso. Ver variable 'data'. Los archivos de entrada (*.txt) tienen que ser 
; matrices de (npts*ntim,4). Las cuatro columnas son, de izquierda a derecha: latitud, longitud, altura (en metros) y 
; emision (en Tn/Yr). 
;-----------------------------------------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;-----------------------------------------------------------------------------------------------------------------------
begin
 npts                                        = 558 ; number of dots
 ntim                                        = 24  ; number of periods (one day).
 nlev                                        = 2  ; number of vertical levels.
 ano                                         = 2023;
 mes                                         = 08;
 dia                                         = 01;
 hora                                        = 00;
 nesp                                        = 30 ; number of chemical species in RADM2 chemical mechanism.
 data                                        = (/"si"   \  ; pm10 available
                                                ,"no"   \  ; co available
                                                ,"no"   \  ; so2 available
                                                ,"no"   \  ; voc available (total)
                                                ,"no"/) \  ; no available
;-----------------------------------------------------------------------------------------------------------------------
 malla                                       = addfile("/data/users/ysilva/WRF_CHEM4/WPS/geo_em.d01.nc","r")
 wrfinput                                    = addfile("/data/users/ysilva/WRF_CHEM4/WRF/run/wrfinput_d01.nc","r")
;-----------------------------------------------------------------------------------------------------------------------
 dx                                          = wrfinput@DX
 dy                                          = wrfinput@DY
;-----------------------------------------------------------------------------------------------------------------------
 zero                                        = new((/npts*ntim,4/),"float")
; [Tn/Yr]->[mol/km2/hr]*10^12/(365*24*M*dx*dy) para gases. M peso molecular.
; [Tn/Yr]->[ug/m2/s]*10^12/(365*24*3600*M*dx*dy) para material particulado.
if(data(0).eq."si")then                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             e_pm10                                     = asciiread("fuegoP_20201019.txt",(/npts*ntim,4/),typeof(zero))
e_pm10(:,3)                                = e_pm10(:,3)*10^9/(24*3600*200*dx*dy)
;  e_pm10(:,3)                             = e_pm10(:,3)*10^9/(dx*dy)  Para SENAMHI
;  e_pm25i                                 = e_pm10*0.67*0.046*0.25
;  e_pm25j                                 = e_pm10*0.67*0.046*0.75
  e_pm25i                                  = e_pm10*0.25
  e_pm25j                                  = e_pm10*0.75
  e_pm_10                                  = e_pm10
  e_eci                                    = e_pm10*0.67*0.029*0.864
  e_ecj                                    = e_pm10*0.67*0.006*0.23
  e_orgi                                   = e_pm10*0.67*0.006*0.77
  e_orgj                                   = e_pm10*0.67*0.40*0.19
  e_so4i                                   = e_pm10*0.67*0.40*0.81
  e_so4j                                   = e_pm10*0.67*0.323*0.94
  e_no3i                                   = e_pm10*0.67*0.323*0.06
  e_no3j                                   = e_pm10*0.33
 else                                             
  e_pm10                                     = zero
  e_pm25i                                    = zero
  e_pm25j                                    = zero 
  e_pm_10                                    = zero
  e_eci                                      = zero
  e_ecj                                      = zero
  e_orgi                                     = zero
  e_orgj                                     = zero
  e_so4i                                     = zero
  e_so4j                                     = zero
  e_no3i                                     = zero
  e_no3j                                     = zero
 end if                                            
;----------------------------------------------------------------------------------------------------------------------
 if(data(4).eq."si")then                          
  e_no                                       = asciiread("no.txt",(/npts*ntim,4/),typeof(zero)) 
  e_no(:,3)                                  = e_no(:,3)*10^9/(dx*dy)
 else                                                
  e_no                                       = zero
 end if                                            
 if(data(1).eq."si")then                              
  e_co                                       = asciiread("co.txt",(/npts*ntim,4/),typeof(zero))
  e_co(:,3)                                  = e_co(:,3)*10^12/(365*24*28*dx*dy)
 else                                               
  e_co                                       = zero
 end if                                             
 if(data(2).eq."si")then                            
  e_so2                                      = asciiread("so2.txt",(/npts*ntim,4/),typeof(zero))
  e_so2(:,3)                                 = e_so2(:,3)*10^9/(dx*dy)
 else                                             
  e_so2                                      = zero
 end if                                       
;----------------------------------------------------------------------------------------------------------------------
 if(data(3).eq."si")then                          
  e_voc                                      = asciiread("voc.txt",(/npts*ntim,4/),typeof(zero))
  e_voc(:,3)                                 = e_voc(:,3)*10^12/(365*24*1*dx*dy)
  e_iso                                      = e_voc
  e_eth                                      = e_voc
  e_hc3                                      = e_voc
  e_hc5                                      = e_voc
  e_hc8                                      = e_voc
  e_xyl                                      = e_voc
  e_ol2                                      = e_voc
  e_olt                                      = e_voc
  e_oli                                      = e_voc
  e_tol                                      = e_voc
  e_csl                                      = zero
  e_hcho                                     = e_voc
  e_ald                                      = e_voc
  e_ket                                      = e_voc
  e_ora2                                     = zero
  e_nh3                                      = e_no*0.034
 else                                             
  e_voc                                      = zero
  e_iso                                      = zero
  e_eth                                      = zero
  e_hc3                                      = zero
  e_hc5                                      = zero
  e_hc8                                      = zero
  e_xyl                                      = zero
  e_ol2                                      = zero
  e_olt                                      = zero
  e_oli                                      = zero
  e_tol                                      = zero
  e_csl                                      = zero
  e_hcho                                     = zero
  e_ald                                      = zero
  e_ket                                      = zero
  e_ora2                                     = zero
  e_nh3                                      = zero
 end if                                           
 ;----------------------------------------------------------------------------------------------------------------------
 dados                                       = (/e_iso,e_so2,e_no,e_co,e_eth,e_hc3,e_hc5,e_hc8,e_xyl,e_ol2       \
                                                ,e_olt,e_oli,e_tol,e_csl,e_hcho,e_ald,e_ket,e_ora2,e_nh3         \
                                                ,e_pm25i,e_pm25j,e_pm_10,e_eci,e_ecj,e_orgi,e_orgj,e_so4i        \
                                                ,e_so4j,e_no3i,e_no3j/)
 espec                                       = (/"E_ISO","E_SO2","E_NO","E_CO","E_ETH","E_HC3","E_HC5","E_HC8"   \
                                                ,"E_XYL","E_OL2","E_OLT","E_OLI","E_TOL","E_CSL","E_HCHO"        \
                                                ,"E_ALD","E_KET","E_ORA2","E_NH3","E_PM25I","E_PM25J","E_PM_10"  \
                                                ,"E_ECI","E_ECJ","E_ORGI","E_ORGJ","E_SO4I","E_SO4J","E_NO3I"    \
                                                ,"E_NO3J"/)
 lat_o                                       = e_pm10(0:npts-1,0)
 lon_o                                       = e_pm10(0:npts-1,1)
 alt_o                                       = e_pm10(0:npts-1,2)
;----------------------------------------------------------------------------------------------------------------------
 XLAT_M                                      = wrf_user_getvar(malla,"XLAT_M",-1)
 XLONG_M                                     = wrf_user_getvar(malla,"XLONG_M",-1)
 dSizes                                      = getfiledimsizes(wrfinput)
 west_east                                   = dSizes(2)
 south_north                                 = dSizes(3)
 wrfchemi                                    = "wrfchemi_d01.nc"
 system("/bin/rm -f "+wrfchemi)                
 f                                           = addfile(wrfchemi,"c")
 setfileoption(f,"DefineMode",True)               
 fAtt                                        = True
 att_wrfinput                                = getvaratts(wrfinput)
 if(.not.all(ismissing(att_wrfinput)))then    
  do i = 0, dimsizes(att_wrfinput)-1             
   if(att_wrfinput(i).ne."TITLE")then           
    fAtt@$att_wrfinput(i)$                   = wrfinput@$att_wrfinput(i)$
   else                                         
    fAtt@$att_wrfinput(i)$                   = "OUTPUT FROM MLV EMISSIONS PREPROCESSOR"
   end if                                         
  end do                                          
 end if                                                 
;-------------------------------------------------------------------------------------------------------------------
 Time                                        = -1                  
 DateStrLen                                  = 19
 emissions_zdim                              = nlev
 Times_string                                = new((/ntim/),"string")
 do i = 0, ntim-1                                 
  if(i.lt.10)then                                 
   Times_string(i)                           = "2023-08-01_0"+i+":00:00"
  else                                             
   Times_string(i)                           = "2023-08-01_"+i+":00:00"
  end if                                           
 end do                                               
 Times_character                             = stringtocharacter(Times_string)
 Times                                       = new((/ntim,DateStrLen/),typeof(Times_character))
 Times(:,:)                                  = Times_character(:,0:DateStrLen-1)
 delete(Times@_FillValue)                              
;-------------------------------------------------------------------------------------------------------------------
 fileattdef(f,fAtt)                                     
 dimNames                                    = (/"Time","DateStrLen","west_east","south_north","emissions_zdim"/)
 dimSizes                                    = (/Time,DateStrLen,west_east,south_north,emissions_zdim/)
 dimUnlim                                    = (/True,False,False,False,False/)
 filedimdef(f,dimNames,dimSizes,dimUnlim)    
 filevardef(f,"Times","character",(/"Time","DateStrLen"/))
;-------------------------------------------------------------------------------------------------------------------
 minlat                                      = min(XLAT_M(0,:,:))
 maxlat                                      = max(XLAT_M(0,:,:)) 
 minlon                                      = min(XLONG_M(0,:,:)) 
 maxlon                                      = max(XLONG_M(0,:,:))
 lat1d                                       = fspan(minlat,maxlat,south_north) ; verificar la cantidad de puntos!
 lon1d                                       = fspan(minlon,maxlon,west_east)
;-------------------------------------------------------------------------------------------------------------------
 hgt                                         = wrfinput->HGT(0,:,:)                    
 ph                                          = wrfinput->PH(0,:,:,:)
 phb                                         = wrfinput->PHB(0,:,:,:)
 geo                                         = (ph+phb)/9.81
 dimgeo                                      = dimsizes(geo)
 z                                           = 0.5*(geo(0:dimgeo(0)-2,:,:)+geo(1:dimgeo(0)-1,:,:))
 emiss                                       = new((/ntim,nlev,south_north,west_east/),typeof(zero))
 emiss!0                                     = "time"
 emiss!1                                     = "lev"
 emiss!2                                     = "lat"
 emiss!3                                     = "lon"  
 emiss&lon                                   = lon1d
 emiss&lat                                   = lat1d
 emiss&lon@units                             = "degrees_east"
 emiss&lat@units                             = "degrees_north"
 taxa                                        = new((/ntim,npts/),typeof(zero))
 lat                                         = new(npts,typeof(zero))
 lon                                         = new(npts,typeof(zero))
 emi                                         = new(npts,typeof(zero))
 do esp = 0, nesp-1                                     
  do tim = 0, ntim-1                               
   taxa(tim,:)                               = dados(esp,npts*tim:npts*(tim+1)-1,3)
   do lev = 0, nlev-1                        
    c                                        = 0
    do npt = 0, npts-1                       
     if(abs(alt_o(npt)-avg(z(lev,:,:)-hgt(:,:))).lt.(0.5*avg(geo(lev+1,:,:)-geo(lev,:,:))))then
      lat(c)                                 = lat_o(npt)
      lon(c)                                 = lon_o(npt)
      emi(c)                                 = taxa(tim,npt)
    else 
      lat(c)                                 = -999
      lon(c)                                 = -999
      emi(c)                                 = -999
     end if                                    
     c                                       = c+1
    end do                                      
    obsij                                    = wrf_user_latlon_to_ij(malla,lat,lon)  
    do pe=0,npts-1                                  
     if(ismissing(emiss(tim,lev,obsij(pe,0),obsij(pe,1)))) then
      emiss(tim,lev,obsij(pe,0),obsij(pe,1)) = emi(pe)
     else                                        
      emiss(tim,lev,obsij(pe,0),obsij(pe,1)) = emiss(tim,lev,obsij(pe,0),obsij(pe,1))+emi(pe)
     end if                                        
    end do                                        
   end do                                           
  end do                                             
  emiss_wrf                                  = rgrid2rcm(emiss&lat,emiss&lon,emiss,XLAT_M(0,:,:),XLONG_M(0,:,:),1)
  emiss_wrf!0                                = "Time"
  emiss_wrf!1                                = "emissions_zdim"
  emiss_wrf!2                                = "south_north"
  emiss_wrf!3                                = "west_east"
  emiss_wrf@_FillValue                       = 0
  delete(emiss_wrf@_FillValue)                  
  emiss_wrf@coordinates                      = "XLONG XLAT"
  emiss_wrf@stagger                          = ""
  if(esp.lt.19)then                              
   emiss_wrf@units                           = "mol km^-2 hr^-1"
  else                                         
   emiss_wrf@units                           = "ug/m3 m/s"
  end if                                             
  emiss_wrf@description                      = "EMISSIONS"
  emiss_wrf@MemoryOrder                      = "XYZ"
  emiss_wrf@FieldType                        = 104
;-------------------------------------------------------------------------------------------------------------------
  filevardef(f,espec(esp),typeof(emiss_wrf),getvardims(emiss_wrf))
  f->$espec(esp)$                            = emiss_wrf
;-------------------------------------------------------------------------------------------------------------------
 end do                                      
 f->Times                                    = Times
 print("**************** Fin ****************")
end
