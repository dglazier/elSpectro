import pdg
api = pdg.connect()

particles = api.get_particles()

for p in particles :
    print(p)


rho = api.get_particles_by_name("rho(770)")
print(rho)
#print(rho[0].mass,rho[0].width)
#print(rho[0].exclusive_branching_fractions )

pi_minus = api.get_particle_by_name('pi-')
#print('MC ID  = ', pi_minus.mcid)
#print('mass   = ', pi_minus.mass, 'GeV')
print('spin J = ', pi_minus.quantum_J)
print('spin P = ', pi_minus.quantum_P)
print('spin C = ', pi_minus.quantum_C)

#import pdg
#api = pdg.connect()
#for bf in api.get_particle_by_name('B+').exclusive_branching_fractions():
#    print('%-60s    %4s    %s' % (bf.description, bf.is_limit, bf.value))

omega = api.get_particle_by_name("omega(782)")
print('MC ID  = ', omega.mcid)
print('mass   = ', omega.mass, 'GeV')
print('spin J = ', omega.quantum_J)
print('spin P = ', omega.quantum_P)
print('spin C = ', omega.quantum_C)
#for bf in omega.exclusive_branching_fractions():
#    print('%-60s    %4s    %s' % (bf.description, bf.is_limit, bf.value))

#for decay in omega.exclusive_branching_fractions():
#    decay_products = [p.item.name for p in decay.decay_products]
#    print(decay_products)
#    if 'pi+' in decay_products:
#        print(format(decay.description,'40s'), decay.display_value_text)

import ROOT
ROOT.gSystem.Load("${ELSPECTRO}/lib/libelSpectro")
ROOT.gROOT.ProcessLine(".x ../ParticleFactory.C")

h_omega = ROOT.phoPro.Hadron("omega",omega.mcid,"Meson")
h_omega.SetMass(omega.mass)
h_omega.SetWidth(omega.width)


h_omega.Print()

#particles.contains("rho")
#print(api.get("omega(782)"))

def exists(pdg_id):
    particles = api.get_particles()
    for p in particles :
        if p.pdgid == pdg_id :
            return True

    return False

print(exists("63623"))
