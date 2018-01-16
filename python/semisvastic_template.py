####################
# configuration file


####################
# general parameters

Nevents = NEVENTS
CutOption = Timing		# Timing->cut every photon PRODUCED after 34 ns , LightColl->cut 99 photons every 100 to speed up the simulation
Seed = SEED

###############
# Crystal
Crystal_x = 11
Crystal_y = 11
Crystal_z = 3
tilt_angle = 0.
surface_type = 5		# Tipe of surface crystal-air (see G4 manual)
				# 0=polished ; 1=polishedfrontpainted ; 2=polishedbackpainted ; 3=ground ; 4=groundfrontpainted ; 5=groundbackpainted 
wrapping_refl = 0.97
SigmaAlpha = 0.
SpecularSpike = 1.
SpecularLobe = 0.
BackScatter = 0.

###############
# SiPM

NSiPM = 1
PDEoption = FBK5x5_20um		# PDE available: KETEK3x3_25um ; KETEK3x3_15um : KETEKDIFFUSE_15um ; FBK5x5_20um ; HPK_MPPC ; PDEFBKDIFFUSE_20um
SiPM_x = 1.
SiPM_y = 1.
SiPM_z = 0.8
SiPM_x_pos=|-5|-4|-3|-2|-1|0|1|2|3|3|3|3|3|3|3|3|5|4|3|2|1|0|-1|-2|-3|-3|-3|-3|-3|-3|-3|-3|
SiPM_y_pos=|-3|-3|-3|-3|-3|-3|-3|-3|-5|-4|-3|-2|-1|0|1|2|3|3|3|3|3|3|3|3|5|4|3|2|1|0|-1|-2|

###############
# Source

SourceDistribution = uniform
Energy = 100
Particle = pi+
X_SourcePosition=|0.|0.|0.|0.|0.|0.|0.|0.|0.|
Y_SourcePosition=|-17.5|-13.4|-8.8|-4.4|0.|+4.4|+8.8|+13.4|+17.5|




