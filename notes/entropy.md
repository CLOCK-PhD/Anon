# ENTROPIE

## Notes de la vidéo "Entropy (for Data Science) Clearyly Explained!!!"

### Introduction

- Très utilisée en data science
- Utilisée pour construire des arbres de classification
- à la base de la *mutual information* qui quantifie la relation entre deux choses
- à la base de l'*entropie relative* et de la *cross entropy" (t-sne et umap)
- L'entropie est utilisée pour **quantifier les similarités et les différences**

### Comment l'entropie quantifie les similarités et les différences

- Pour parler de l'entropie, il faut d'abord comprendre le principe de *surprise*

####  Surprise

Exemple des boules noires et blanches dans un sac:
- On a trois sacs avec des boules noires et blanches
- Le premier a beaucoup plus de blanches que de noires
- Le deuxième a beaucoup plus de noires que de blanches
- Le troisième a autant de blanches que de noires

- Piocher une boule blanche dans le premier sac est le résultat auquel on s'attend le plus
- Piocher une boule noire dans le premier sac est un événement surprenant puisqu'on a beaucoup moins de chances de la piocher
- C'est l'inverse pour le deuxième sac, où il est plus facile de piocher une boule noire qu'une boule blanche
- Dans le dernier sac, on a autant de chance de piocher une boule noire ou une boule blanche : la surprise est la égale dans les deux cas

La suprise est en quelque sorte inversement liée à la probabilité.
En conséquence, quand la probabilité de piocher une boule est faible,la surprise est élevée

##### Calcul de la surprise

Ce n'est pas l'inverse de la probabilité (1/P).

Avec l'exemple d'un lancé de pièce: Imaginons qu'on lance plusieurs fois une pièce, et qu'elle donne systématiquement "face".
Si on fait un lancer supplémentaire, serait-on surpris si on avait encore "face" ? Pas du tout, donc on pourrait estimer la surprise à 0, alors que $P(face) = 1$. Donc $1/P(face) = 1/1 = 1$.

C'est pourquoi, pour calculer la surprise, on utilise le log de l'inverse de la probabilité : $\log(1/p(face))$
Et là : $\log(1/p(face)) = \log(1/1) = 0$, la surprise est donc de 0.

À l'opposé, puisque la probabilité d'obtenir pile est de 0, et qu'on n'en aura donc jamais, cela n'a pas de sens de quantifier la surprise d'un événement qui ne peut jamais se produire.

On a donc : $\log(1/p(pile)) = \log(1/0) = \log(1) - \log(0) =$ non défini, car log(0) n'est pas défini. Et c'est OK parce qu'on parle de la surprise pour un événement qui ne peut jamais se produire.

Si on a deux outputs (pile et face), on utilise un log base 2 pour les calculs

Maintenant qu'on a défini ce qu'était la surprise, imaginons que notre pièce donne Face 90% du temps et pile 10% du temps.

Calculons la surprise :
- Surprise pile : $log_2(\frac{1}{p(face)}) = log_2(\frac{1}{0.9}) = log_2(1) - log2(0.9) = 0.15$
- Surprise face : $log_2(\frac{1}{p(pile)}) = log_2(\frac{1}{0.1}) = log_2(1) - log_2(0.1) = 3.32$

Puisqu'obtenir face est plus rare qu'obtenir pile, la surprise pour pile est beaucoup plus élevée.

Maintenant, faisons 3 lancers  successifs : face, face, pile
La proba d'obtenir 2 faces et 1 pile est $0.9  \times  0.9  \times  031$

Si on veut estimer la surprise, on place cette probabilité dans l'équation de la surprise :
$Surpise = log_2 (1/0.9 \times 0.9 \times 0.1) = \log_2(1) - log_2(0.9  \times  0.9  \times  0.1) = \log_2(1) - [log_2(0.9)+ log_2(0.9)+ log_2(0.1)] = 0 - 0.15 + 0.15 + 0.32 = 3.62$

Le point important ici est qu'on observe que la Surprise totale pour un lancer de pièces est la **somme des surprises de chaque lancer individuel**

##### De la surprise à l'entropie
On a donc :

                face    pile

    p(x)        0.9     0.1

    Surprise    0.15    3.32

Maintenant, si on veut estimer la surprise totale après 100 lancers successifs, on va faire une approximation du nombre de fois qu'on va avoir face en multiplant la probabilité par le nombre de lancers : $0.9 * 100$, et on estime la surprise d'obtenir face en multipliant par 0.15.

Ainsi, $(0.9  \times  100)  \times  0.15$ représente la surprise attendue d'avoir face dans 100 lancers.

De la même manière, on peut approcher le nombre attendu de pile en multipliant la proba de pile par le nombre de lancers : $0.1  \times  100$, et on estime la surprise d'obtenir pile en multipliant par 3.32.

Ainsi, $(0.1  \times  100)  \times  3.32$ représente la surprise attendue d'avoir pile dans 100 lancers.

Pour trouver la surprise totale, on additionne les deux termes : 

$(0.9  \times  100)  \times  0.15 + (0.1  \times  100)  \times  3.32 = 46.7$

On vient donc d'estimer la surprise totale pour 100 lancers successifs, mais qu'en est-il de l'entropie ?

Si on divise tout par le nombre de lancers (100), on obtient alors la surprise moyenne par lancer : $46.7 / 100 = 0.47$

Ainsi, en moyenne, on s'attend à ce que la surprise soit de 0.47 à chaque fois qu'on lance la pièce. **Ce qui nous donne l'entropie de la pièce : la surprise attendue à chaque fois qu'on lance la pièce**.

En statistiques, on dit que l'entropie est la *valeur attendue* de la surprise.

$E(surprise) = (0.9 * 100) * 0.15 + (0.1 * 100) * 3.32 / 100$

On peut simplifier en retirant les 100:
$E(surprise) = (0.9* 0.15) + (0.1 * 3.32) = 0.47$

$$E(surprise) = \sum x P(X=x)$$
Avec $x$ la valeur spécifique de la surprise, et  $P(X=x)$ la probabilité d'observer cette valeur spécifique pour la surprise.

L'entropie est donc *la surprise moyenne à laquelle on peut s'attendre*, il est possible de la dériver à partir de la surprise.

En effet, si on implante l'équation de la surprise pour la valeur $x$ :
$$\sum \log \frac{1}{p(x)}P(X=x)$$

Puis qu'on implante la probabilité, on retrouve la formule de l'entropie :

$$H = \sum \log \frac{1}{p(x)}p(x)$$

avec $\log \frac{1}{p(x)}$ la surprise et $p(x)$ la probabilité de la surprise.

Cette equation n'est pas la forme habituelle de l'entropie qu'on peut trouver.

En changeant l'ordre des termes :

$$H = \sum p(x)\log \frac{1}{p(x)}$$

Puis en utilisant les propriétés des logs pour convertir la fraction en soustraction :

$$H = \sum p(x) [\log(1) - \log(p(x))]$$

Comme $\log(1) = 0$, on a :

$$H = \sum p(x) [\log(1) - \log(p(x))]$$

En multipliant les deux termes par $p(x)$:

$$H = \sum -p(x)\log(p(x))$$

En retirant le signe moins de la somme, on se retrouve enfin avec la formule de l'entropie telle que décrite par Claude Shannon en 1948 :

$$H = - \sum p(x)\log(p(x))$$

### Calcul de l'entropie des boules noires et blanches

En retournant sur notre premier exemple avec les sacs de boules.

Si notre premier sac contient 6 boules blanches et une noire, on aura pour les probabilités:
- p(boule blanche) = 6/7
- p(boule noire) = 1/7

l'entropie sera
- pour les boules blanches : $6/7 \times \log_2\frac{1}{6/7}$
- pour les boules noires : $1/7  \times  log_2\frac{1}{1/7}$

L'entropie est donc : $H = (0.86  \times  0.22) + (0.14  \times  2.81) = 0.59$

Même si la surprise associée au tirage d'une boule blanche (0.22) est plus petite que celle associée à celui d'une boule noire (2.81), la probabilité de tirer une boule blanche (0.86) est beaucoup plus élevée que de tirer une boule noire (0.14). Ainsi, l'entropie totale (0.59) est beaucoup plus proche de la surprise associée aux boules blanches (0.22) qu'aux boules noires (2.81).

Pour le deuxième sac qui contient plus de boules noires que de blanches, disons 10 boules noires et une boule blanche, on aura donc comme probabilités pour le tirage :
- p(boule blanche) = 1/11
- p(boule noire) = 10/11

L'entropie est donc :
- H(boule blanche) = $1/11  \times  log_2\frac{1}{1/11}$
- H(boule noire) = $10/11  \times  log_2\frac{1}{10/11}$

L'entropie totale est : $(0.09 \times 3.46) + (0.91 \times 0.14) = 0.44$

Dans ce cas, la surprise de tirer une boule blanche est relativement élevée (3.46), mais la probabilité que cela se produise est tellement basse que l'entropie totale (0.44) est beaucoup plus proche de la surprise associée au tirage d'une boule noire (0.14).

On peut également voir que la valeur de l'entropie (et donc la surprise attendue), est plus faible pour le sac 2 que pour le sac 1. Cela a du sens dans la mesure où dans le sac 2, on a une plus grand chance de tirer une boule avec une surprise plus faible.

Enfin, considérons le dernier sac, où on a un nombre égal de boules noires et de boules blanches, disons 7 pour chaque cas.

On a ici un cas particulier, puisque les probabilités de tirer une boule noire ou blanche sont égales $p(boule blanche) = 7/14 = p(boule noire)$.

Ainsi, le calcul de l'entropie donne :

$H = 7/14  \times  log_2(\frac{1}{7/14}) + 7/14  \times  log_2(\frac{1}{7/14}) = (0.5 \times 1) + (0.5 \times 1)$

Et donc : $H = 1$

Dans ce cas, on retrouve la valeur maximale d'entropie qu'il est possible d'avoir. Dans ce cas, même si la surprise de tirer une boule noire ou une boule blanche est modérée (1), on obtient toujours la même surprise à chaque tirage qui n'est jamais pondérée par des valeurs inférieures de surprise comme dans les deux cas précédents.

Ainsi, on peut utiliser l'entropie  pour quantifier la similarité ou la différence dans le nombre de boules noires ou blanches dans chaque sac.

L'entropie est la plus haute quand on a le même nombre de chaque type de boule, et elle diminue avec l'augmentation de la différence entre le nombre de boules de chaque type.

L'entropie est décrite par la courbe suivante : voir Shannon p11.