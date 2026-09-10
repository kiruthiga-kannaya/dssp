# Création du dictionnaire contenant les coordonnées C O N H de chaque residu

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

##############################
######## HELICE ALPHA ########
##############################

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


# Pour ne garder que les segments de deux turns minimum
helice_alpha = []

for segment in segment_4_turn:
    if segment["fin"] != segment["debut"] + 4 :
        helice_alpha.append(segment)

print(helice_alpha)

##############################
######## HELICE 3 10 ########
##############################

# Repérage du motif 3-turn
liaison_3_turn = []

for liaison in liaison_h :
    if liaison["donneur"] == liaison["accepteur"] + 3 :
        liaison_3_turn.append(liaison)


print(liaison_3_turn)


# Identification des hélices 3_10

segment_3_turn = []

if len(liaison_3_turn) > 0:

    debut = liaison_3_turn[0]["accepteur"]

    for i in range(1, len(liaison_3_turn)) : 

        actuel = liaison_3_turn[i]["accepteur"]
        precedent = liaison_3_turn[i-1]["accepteur"]

        if actuel - precedent != 1 :

            # alors l'enchainement précédent se termine avec le donneur de la liaison précédente
            segment_3_turn.append({"debut" : debut, "fin" : liaison_3_turn[i-1]["donneur"]})

            # enchainement de l'hélice prochain
            debut = actuel
    
    # debut de l'hélice et la fin est le donneur du dernier aa de la liste de dictionnaire
    segment_3_turn.append({"debut" : debut, "fin" : liaison_3_turn[-1]["donneur"]})


# Pour ne garder que les segments de deux turns minimum
helice_3_10 = []

for segment in segment_3_turn:
    if segment["fin"] != segment["debut"] + 3 :
        helice_3_10.append(segment)

print(helice_3_10)

###############################
######## HELICE 5 (pi) ########
###############################

# Repérage du motif 5-turn
liaison_5_turn = []

for liaison in liaison_h :
    if liaison["donneur"] == liaison["accepteur"] + 5 :
        liaison_5_turn.append(liaison)


print(liaison_5_turn)


# Identification des hélices pi

segment_5_turn = []

if len(liaison_5_turn) > 0:

    debut = liaison_5_turn[0]["accepteur"]

    for i in range(1, len(liaison_5_turn)) : 

        actuel = liaison_5_turn[i]["accepteur"]
        precedent = liaison_5_turn[i-1]["accepteur"]

        if actuel - precedent != 1 :

            # alors l'enchainement précédent se termine avec le donneur de la liaison précédente
            segment_5_turn.append({"debut" : debut, "fin" : liaison_5_turn[i-1]["donneur"]})

            # enchainement de l'hélice prochain
            debut = actuel
    
    # debut de l'hélice et la fin est le donneur du dernier aa de la liste de dictionnaire
    segment_5_turn.append({"debut" : debut, "fin" : liaison_5_turn[-1]["donneur"]})


# Pour ne garder que les segments de deux turns minimum
helice_pi = []

for segment in segment_5_turn:
    if segment["fin"] != segment["debut"] + 5 :
        helice_pi.append(segment)

print(helice_pi)


###############################
######## FEUILLET BÊTA ########
###############################

# Bridge (pont) --> ladder (echelle) --> sheet (feuillet)

#################
#### BRIDGES ####
#################

## PARALELLES BRIDGES ##
########################

# fonction pour savoir s'il existe une liaison hydrogène entre un residu 1 et 2 (permet de faciliter le code)
def is_hbond(residu_1, residu_2):

    for liaison in liaison_h : 
        if liaison["accepteur"] == residu_1 and liaison["donneur"] == residu_2 : 
            return True 

    else :
        return False


# sequence avec le numéro de chaque résidus
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
###############################

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


#################
#### LADDERS ####
#################

### PARALELLES LADDERS ###
##########################

ladder_parallele = []

lettre_ladder = "a"

if len(bridge_parallele) > 0:

    debut_1 = bridge_parallele[0]["residu_1"]
    debut_2 = bridge_parallele[0]["residu_2"]

    for i in range(1, len(bridge_parallele)) : 

        actuel_1 = bridge_parallele[i]["residu_1"]
        precedent_1 = bridge_parallele[i-1]["residu_1"]

        actuel_2 = bridge_parallele[i]["residu_2"]
        precedent_2 = bridge_parallele[i-1]["residu_2"]

        if actuel_1 -precedent_1 == 1 and actuel_2 - precedent_2 == 1:
            continue
        else : 
            #on met fin à la ladder
            ladder_parallele.append({lettre_ladder : {"sequence_1" : [debut_1, precedent_1], "sequence_2" : [debut_2, precedent_2]} }) 
            
            #nouveau ladder
            debut_1 = actuel_1
            debut_2 = actuel_2

            #lettre ladder suivant
            lettre_ladder = chr(ord(lettre_ladder) + 1)

    # Ajouter la dernière ladder
    ladder_parallele.append({lettre_ladder : {"sequence_1": [debut_1, bridge_parallele[-1]["residu_1"]],
                             "sequence_2": [debut_2, bridge_parallele[-1]["residu_2"]]}})

print(ladder_parallele)

### ANTI-PARALELLES LADDERS ###
##############################

ladder_anti_parallele = []

lettre_ladder = "A"

if len(bridge_anti_parallele) > 0:

    debut_1 = bridge_anti_parallele[0]["residu_1"]
    debut_2 = bridge_anti_parallele[0]["residu_2"]

    for i in range(1, len(bridge_anti_parallele)) : 

        actuel_1 = bridge_anti_parallele[i]["residu_1"]
        precedent_1 = bridge_anti_parallele[i-1]["residu_1"]

        actuel_2 = bridge_anti_parallele[i]["residu_2"]
        precedent_2 = bridge_anti_parallele[i-1]["residu_2"]

        if actuel_1 -precedent_1 == 1 and actuel_2 - precedent_2 == -1:
            continue
        else : 
            #on met fin à la ladder
            ladder_anti_parallele.append({lettre_ladder : {"sequence_1" : [debut_1, precedent_1], "sequence_2" : [debut_2, precedent_2]} }) 
            
            #nouveau ladder
            debut_1 = actuel_1
            debut_2 = actuel_2

            #lettre ladder suivant
            lettre_ladder = chr(ord(lettre_ladder) + 1)

    # Ajouter la dernière ladder
    ladder_anti_parallele.append({lettre_ladder : {"sequence_1": [debut_1, bridge_anti_parallele[-1]["residu_1"]],
                             "sequence_2": [debut_2, bridge_anti_parallele[-1]["residu_2"]]}})

print(ladder_anti_parallele)


#################
##### SHEET #####
#################

# sheet = ensemble d'une ou plusieurs ladders connectées par des rédidus communs
# un résidu commun : 
### entre 2 ladders au minimum (feuillet beta minimal) 
### 3 ladders (grand feuillet)


# rassembler les ladders parallèles et anti-paralleles 

all_ladders = []

for ladder in ladder_parallele:
    all_ladders.append(ladder)

for ladder in ladder_anti_parallele:
    all_ladders.append(ladder)


# je compare chaque ladders deux à deux et je regarde s'il y a au moins un residu commun entre les deux ladders
# je les regroupe ensemble 

# CONNEXIONS ENTRE LADDERS #
############################

connexions = []

for i in range(len(all_ladders)):

    # recuperer la lettre de la premiere ladder
    lettre_1 = list(all_ladders[i].keys())[0]

    sequence_1 = all_ladders[i][lettre_1]["sequence_1"]
    sequence_2 = all_ladders[i][lettre_1]["sequence_2"]

    # transformer les intervalles de la première ladder en liste de résidus
    residus_ladder_1 = []

    for residu in range(sequence_1[0], sequence_1[1] + 1):
        residus_ladder_1.append(residu)

    if sequence_2[0] <= sequence_2[1]:

        for residu in range(sequence_2[0], sequence_2[1] + 1):
            residus_ladder_1.append(residu)

    else:

        for residu in range(sequence_2[0], sequence_2[1] - 1, -1):
            residus_ladder_1.append(residu)


    # comparer avec les ladders suivantes
    for j in range(i + 1, len(all_ladders)):

        # récupérer la lettre de la deuxième ladder
        lettre_2 = list(all_ladders[j].keys())[0]

        sequence_3 = all_ladders[j][lettre_2]["sequence_1"]
        sequence_4 = all_ladders[j][lettre_2]["sequence_2"]

        # transformer les intervalles de la deuxième ladder en liste de résidus
        residus_ladder_2 = []

        for residu in range(sequence_3[0], sequence_3[1] + 1):
            residus_ladder_2.append(residu)

        if sequence_4[0] <= sequence_4[1]:

            for residu in range(sequence_4[0], sequence_4[1] + 1):
                residus_ladder_2.append(residu)

        else:

            for residu in range(sequence_4[0], sequence_4[1] - 1, -1):
                residus_ladder_2.append(residu)


        # chercher les résidus communs
        residus_communs = []

        for residu in residus_ladder_1:

            if residu in residus_ladder_2:

                if residu not in residus_communs:
                    residus_communs.append(residu)


        # si les deux ladders ont au moins un résidu commun
        if len(residus_communs) > 0:

            connexions.append({
                "ladder_1": lettre_1,
                "ladder_2": lettre_2,
                "residus_communs": residus_communs
            })


print(connexions)




# Rassembler les connextions #
##############################

Sheets = {}

nom_sheet = "A"

for connexion in connexions:

    ladder1 = connexion["ladder_1"]
    ladder2 = connexion["ladder_2"]

    trouve_1 = False
    trouve_2 = False

    sheet_1 = ""
    sheet_2 = ""

    # Chercher dans quels sheets se trouvent les deux ladders
    for sheet in Sheets:

        if ladder1 in Sheets[sheet]:
            trouve_1 = True
            sheet_1 = sheet

        if ladder2 in Sheets[sheet]:
            trouve_2 = True
            sheet_2 = sheet

    # Cas 1 aucune des deux ladders n'est encore dans un sheet
    if not trouve_1 and not trouve_2:

        Sheets[nom_sheet] = [ladder1, ladder2]

        nom_sheet = chr(ord(nom_sheet) + 1)

    # Cas 2 ladder1 est deja dans un sheet
    elif trouve_1 and not trouve_2:

        Sheets[sheet_1].append(ladder2)

    # Cas 3 ladder2 est deja dans un sheet
    elif not trouve_1 and trouve_2:

        Sheets[sheet_2].append(ladder1)

    # Cas 4 les deux ladders sont deja dans le meme sheet
    elif sheet_1 == sheet_2:

        pass

    # Cas 5 les deux ladders sont dans deux sheets differents
    else:

        for ladder in Sheets[sheet_2]:

            if ladder not in Sheets[sheet_1]:
                Sheets[sheet_1].append(ladder)

        del Sheets[sheet_2]


print(Sheets)



#################
#### SUMMARY ####
#################
summary = {}

# Tous les résidus commencent sans structure
for residu in coordonnee_proteine:
    summary[residu] = ""

##################### helices

for helice in helice_alpha:
    debut = helice["debut"]
    fin = helice["fin"]

    for residu in range(debut, fin + 1):
        summary[residu] = "H"


##################### beta-sheets

for ladder in all_ladders:

    for lettre_ladder in ladder:

        sequence_1 = ladder[lettre_ladder]["sequence_1"]
        sequence_2 = ladder[lettre_ladder]["sequence_2"]

        debut_1 = sequence_1[0]
        fin_1 = sequence_1[1]

        debut_2 = sequence_2[0]
        fin_2 = sequence_2[1]

        # Longueur de la ladder
        longueur = abs(fin_1 - debut_1) + 1

        # Ladder d'un seul bridge
        if longueur == 1:
            structure = "B"

        # Ladder de plusieurs bridges
        else:
            structure = "E"

        # 1ere sequence
        for residu in range(debut_1, fin_1 + 1):
            summary[residu] = structure

        # 2e sequence
        if debut_2 <= fin_2:

            for residu in range(debut_2, fin_2 + 1):
                summary[residu] = structure

        else:

            for residu in range(debut_2, fin_2 - 1, -1):
                summary[residu] = structure


## TABLEAU FINAL 

import pandas as pd

# Numeros des residus
numero_residu = list(summary.keys())

# Noms des residus
nom_residu = []

for residu in numero_residu:
    nom_residu.append(coordonnee_proteine[residu]["nom"])

# Structures secondaires
structure = list(summary.values())

# dataframe
df = pd.DataFrame()

df["numero_residu"] = numero_residu
df["nom_residu"] = nom_residu
df["structure"] = structure

print(df)

df_transpose = df.T
print(df_transpose)

df_transpose = df.T

# pour que ce soit plus visuel
taille = 20

for debut in range(0, len(df), taille):

    fin = debut + taille

    morceau = df[debut:fin]

    print(morceau.T)
    print()

#enregistrer le tableau dans le dossier results
df.to_csv("results/summary.csv", index=False)