# Extraire le nombre d'acide aminé dans la protéine, pour structurer le dictionnaire : numéro de résidu + nom du résidu
coordonnee_proteine = {}

with open("data/1P0R.pdb", "r") as proteine:
    for ligne in proteine:

        if ligne.startswith("ATOM"):

            numero_residu = int(ligne[22:26])
            nom_residu = ligne[17:20].strip()
            nom_atome = ligne[12:16].strip()

            if numero_residu not in coordonnee_proteine:
                coordonnee_proteine[numero_residu] = {}
                coordonnee_proteine[numero_residu]["nom"] = nom_residu

# extraire les coordonnées de C, O, N, H

            if nom_atome == "N":
                coordonnee_proteine[numero_residu]["N"] = [
                    float(ligne[30:38]),
                    float(ligne[38:46]),
                    float(ligne[46:54])
                ]
            
            if nom_atome =="C":
                coordonnee_proteine[numero_residu]["C"] = [
                    float(ligne[30:38]),
                    float(ligne[38:46]),
                    float(ligne[46:54])
                ]
            
            if nom_atome =="O":
                coordonnee_proteine[numero_residu]["O"] = [
                    float(ligne[30:38]),
                    float(ligne[38:46]),
                    float(ligne[46:54])
                ]
            
            if nom_atome =="H":
                coordonnee_proteine[numero_residu]["H"] = [
                    float(ligne[30:38]),
                    float(ligne[38:46]),
                    float(ligne[46:54])
                ]

print(coordonnee_proteine)

print(len(coordonnee_proteine))


def calcule_hbond(coord_accepteur, coord_donneur) : 
    import math
    
     
# Calcul des distances
    xc = coord_accepteur["C"][0]
    yc = coord_accepteur["C"][1]
    zc = coord_accepteur["C"][2]

    xo = coord_accepteur["O"][0]
    yo = coord_accepteur["O"][1]
    zo = coord_accepteur["O"][2]


    xn = coord_donneur["N"][0]
    yn = coord_donneur["N"][1]
    zn = coord_donneur["N"][2]

    xh = coord_donneur["H"][0]
    yh =coord_donneur["H"][1]
    zh = coord_donneur["H"][2]

# les distances sont en angstrom 
    r_ON = math.sqrt((xo - xn)**2 + (yo - yn)**2 + (zo - zn)**2)
    r_OH = math.sqrt((xo - xh)**2 + (yo - yh)**2 + (zo - zh)**2)

    r_CN = math.sqrt((xc - xn)**2 + (yc - yn)**2 + (zc - zn)**2)
    r_CH = math.sqrt((xc - xh)**2 + (yc - yh)**2 + (zc - zh)**2)

    q1 = 0.42
    q2 = 0.20

# calcul de l'énergie
    energie = q1 * q2 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)*332 # en kcal/mole

    return energie 


# print(calcule_hbond(coordonnee_proteine[2], coordonnee_proteine[3]))

# Calcul des liaisons hydrogènes pour chaque paires de résidus  
liaison_h = []

seuil = -0.5

for residu_accepteur in coordonnee_proteine :  

    for residu_donneur in coordonnee_proteine : # pour qu'un residu ne soit pas comparé à lui-même (1-1, 2-2 ...)
        if residu_accepteur == residu_donneur :
            continue 

        if "C" in coordonnee_proteine[residu_accepteur] and "O" in coordonnee_proteine[residu_accepteur] and "N" in coordonnee_proteine[residu_donneur] and "H" in coordonnee_proteine[residu_donneur] :
            resultat = calcule_hbond(coordonnee_proteine[residu_accepteur], coordonnee_proteine[residu_donneur])
    
            if resultat < seuil :
                liaison_h.append({"accepteur" : residu_accepteur, "donneur" : residu_donneur, "energie" :resultat})
        
print(liaison_h)

# Repérage du motif 4-turn
liaison_4_turn = []

for liaison in liaison_h :
    if liaison["donneur"] == liaison["accepteur"] + 4 :
        liaison_4_turn.append(liaison)

print(liaison_4_turn)

# Identification des hélices alphas
helice_alpha = []

for liaison in liaison_4_turn:
    enchainement_aa = []
    enchainement_aa.append(liaison["accepteur"])

    aa_helice : []
    for numero in enchainement_aa:
        if numero+1 - numero == 1 :
            aa_helice[numero, numero +1]

print(aa_helice)
print(helice_alpha)


