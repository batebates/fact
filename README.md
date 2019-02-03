**Outil de Factorisation**

Cet outil permet à l'utilisateur de factoriser une liste de nombre contenu dans un fichier.

**Installation**

Télécharger l'outil et lancer la commande `make fact`


**Utilisation**


Usage: fact -f [FILE] [OPTION]

  Par défault Factorisation optimisé en fonction de la grandeur du nombre
  -d Factorisation par divisions successives
  -p Factorisation par la méthode p-1 de Pollard 
  -r Factorisation par la méthode rho de Pollard 
  -b Factorisation par la méthode rho de Pollard variante de Brent

Exemple:  
  `./fact -f 10chiffres.txt`
  `./fact -f 70chiffres.txt -p`
