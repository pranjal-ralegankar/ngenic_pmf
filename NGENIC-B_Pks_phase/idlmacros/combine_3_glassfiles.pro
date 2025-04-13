
F1 = "../glass1/snapshot_000"
F2 = "../glass2/snapshot_000"
F3 = "../glass3/snapshot_000"

Fout = "combined_glasses.dat"

mfrac1 = 0.15d
mfrac2 = 0.425d
mfrac3 = 0.425d


UnitLength_in_cm =        3.085678d21        ;  1.0 kpc 
UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec 
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
GRAVITY   = 6.672d-8
HUBBLE =  3.2407789d-18   

G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
H0 = HUBBLE * UnitTime_in_s



npart=lonarr(6)
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)
flag_cooling= 0L
num_files= 1L
BoxSize = 0.0D

bytesleft=120
la=intarr(bytesleft/2)


openr,1, F1,/f77_unformatted
readu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
N1=  npart(1)
pos1 = fltarr(3, N1)
readu,1, pos1
id1 = lindgen(N1)+1
close,1

openr,1, F2,/f77_unformatted
readu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
N2=  npart(1)
pos2 = fltarr(3, N2)
id2 = lindgen(N2)+1+N1
readu,1, pos2
close,1

openr,1, F3,/f77_unformatted
readu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
N3=  npart(1)
pos3 = fltarr(3, N3)
id3 = lindgen(N3)+1+N1+N2
readu,1, pos3
close,1


npart(*)=0
npartall(*)=0

npart(1) = N1
npart(2) = N2
npart(3) = N3

npartall(1) = N1
npartall(2) = N2
npartall(3) = N3



M1 = mfrac1 * 3 * H0^2/(8*!DPI*G) * Boxsize^3 / npart(1)
M2 = mfrac2 * 3 * H0^2/(8*!DPI*G) * Boxsize^3 / npart(1)
M3 = mfrac3 * 3 * H0^2/(8*!DPI*G) * Boxsize^3 / npart(1)


massarr(*)=0
massarr(1) = M1
massarr(2) = M2
massarr(3) = M3



vel = fltarr(3,N1+N2+N3)

openw,1, Fout,/f77_unformatted
writeu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
writeu,1, pos1,pos2,pos3
writeu,1, vel
writeu,1, id1,id2,id3
close,1






ind=where(pos1(2,*) lt boxsize/64)
plot,  pos1(0,ind), pos1(1,ind), psym=3

ind=where(pos2(2,*) lt boxsize/64)
oplot,  pos2(0,ind), pos2(1,ind), psym=3, color=255


end
