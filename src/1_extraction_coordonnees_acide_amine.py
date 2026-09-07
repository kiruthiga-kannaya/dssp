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

################
# HELICE ALPHA #
################

# Repérage du motif 4-turn
liaison_4_turn = []

for liaison in liaison_h :
    if liaison["donneur"] == liaison["accepteur"] + 4 :
        liaison_4_turn.append(liaison)


print(liaison_4_turn)


# Identification des hélices alphas

segment_4_turn = []

if len(liaison_4_turn) > 0:

    debut = liaison_4_turn[0]["accepteur"]

    for i in range(1, len(liaison_4_turn)) : 

        actuel = liaison_4_turn[i]["accepteur"]
        precedent = liaison_4_turn[i-1]["accepteur"]

        if actuel - precedent != 1 :

            # alors l'enchainement précédent se termine avec le donneur de la liaison précédente
            segment_4_turn.append({"debut" : debut, "fin" : liaison_4_turn[i-1]["donneur"]})

            # enchainement de l'hélice prochain
            debut = actuel
    
    # debut de l'hélice et la fin est le donneur du dernier aa de la liste de dictionnaire
    segment_4_turn.append({"debut" : debut, "fin" : liaison_4_turn[-1]["donneur"]})


# pour ne garder que les segments de deux turns minimum
helice_alpha = []

for segment in segment_4_turn:
    if segment["fin"] != segment["debut"] + 4 :
        helice_alpha.append(segment)

print(helice_alpha)


#################
# FEUILLET BÊTA #
#################

# coordonnées_protéine, nom et numero residu
# deux boucles pour parcourir 

#pont 
# trio, et que entre deux trio pas de chevauchements -> condition 1
# calcul des liaisons hydrogènes 


# Bridge (pont) --> ladder (echelle) --> sheet (feuillet)

#### BRIDGES ####

### PARALELLES BRIDGES ###

# prendre une fenetre de 3 pour accepteur et une fenetre de 3 pour donneur de la liste de dictionnaire liaison_h
# il faut pas qu'ils se chevauchent 
# voir s'il existe 2 Hbond entre ces fenetres à chaque fois (ca et ca OU ca et ca)
# parcourir les numeros de résidus dans coordonnées protéines 
# je prend numero residu i et j 
# regarde dans liaison_h s'il existe une relation entre accepteur i(formule) et donneur j(formule)
# si il existe une liaison_h entre accepteur i(formule) et donneur j(formule) alors  ajouter la liaison dans bridge parallèles

def is_hbond(residu_1, residu_2):
    for liaison in liaison_h : 
        if liaison["accepteur"] == residu_1 and liaison["donneur"] == residu_2 : 
            return True 
    
    else :
        return False


#sequence avec le numéro de chaque résidus
seq_numero_residu = [] 

for numero_residu in coordonnee_proteine :
    seq_numero_residu.append(numero_residu)


bridge_parallele = []
for i in range(1, len(seq_numero_residu)- 1): 

    for j in range(4, len(seq_numero_residu) -1):

        if i > j : # pour eviter d'avoir des doublons de paires
            continue

        if abs(i - j) < 3: # pour eviter les chevauchements
            continue

        if is_hbond (i- 1, j) and is_hbond (j, i+1) or is_hbond (j-1, i) and is_hbond (i, j+1) :
            bridge_parallele.append({"residu_1" : i, "residu_2" : j})

print(bridge_parallele)


### ANTI-PARALELLES BRIDGES ###

bridge_anti_parallele = []
for i in range(1, len(seq_numero_residu)- 1): 

    for j in range(4, len(seq_numero_residu) -1):

        if i > j : # pour eviter d'avoir des doublons de paires
            continue
        
        if abs(i - j) < 3: # pour eviter les chevauchements
            continue

        if is_hbond (i, j) and is_hbond (j, i) or is_hbond (i-1, j+1) and is_hbond (j-1, i+1) :
            bridge_anti_parallele.append({"residu_1" : i, "residu_2" : j})

print(bridge_anti_parallele)

