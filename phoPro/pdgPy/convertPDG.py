import pdg
api = pdg.connect()

import ROOT
ROOT.gSystem.Load("${ELSPECTRO}/lib/libelSpectro")
ROOT.gROOT.ProcessLine(".x ../ParticleFactory.C+")
outFile = ROOT.TFile("elspectroPDG.root","recreate")

def exists(pdg_id):
    particles = api.get_particles()
    for p in particles :
        if p.pdgid == pdg_id :
            return True

    return False

def is_pipi(decay):
    products = decay.decay_products
    #print(dir(products))
    if len(products) < 2 :
        return False
    if products[0] and products[1] :
        if products[0].item.name=='pi' and products[1].item.name=='pi' :
            return True
    return False

#particles = ['omega(782)','phi(1020)']
particles = ['rho(770)0']
#particles = ['a_0(980)+']
#particles = ['f_2(1270)']
particles = ['a_1(1260)+']
for partname in particles :
    part = api.get_particle_by_name(partname)
    print(part)
    
    h_part = ROOT.phoPro.Hadron(part.mcid,"Meson")
    h_part.SetMass(part.mass)
    if part.width :
        h_part.SetWidth(part.width)
    
    for decay in part.exclusive_branching_fractions():
        print(decay)
        prod_names = []
        prod_mcpids = []
        if is_pipi(decay) :
            prod_names.append(['pi+','pi-'])
            prod_names.append(['pi0','pi0'])
            prod_mcpids.append([211,-211])
            prod_mcpids.append([111,111])
            
        for idec,dec in enumerate(prod_names):
            n_ambigs=len(prod_names)
            decay_hadrons = []
            print('idec',idec)
            for iprod,prod in enumerate(prod_names[idec]):
                print('name,id',prod_names[idec][iprod],prod_mcpids[idec][iprod])
                decay_hadrons.append(ROOT.phoPro.Hadron(prod_names[idec][iprod],prod_mcpids[idec][iprod]))
                
            h_part.AddDecay(decay.value/n_ambigs,decay_hadrons)

            
 
    h_part.Print()
    h_part.Write()
   

    
