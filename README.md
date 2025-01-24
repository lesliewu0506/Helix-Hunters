# Protein Pow(d)er
## Introductie case
Eiwitten zijn essentiële moleculen in het menselijk lichaam, opgebouwd uit lange ketens van aminozuren. De specifieke vouwing van deze eiwitten bepaalt hun functie. Verkeerd gevouwen eiwitten kunnen leiden tot ernstige ziekten zoals kanker, Alzheimer en taaislijmziekte. Het doel van dit project is om algoritmen te ontwerpen en te implementeren die eiwitstructuren optimaliseren door de meest stabiele vouwing te vinden.

### Vouwing
Dit project maakt gebruik van een 3D-grid (optioneel ook 2D) om eiwitten te modelleren.. Hierbij wordt elke aminozuur stap voor stap geplaatst in de grid. Bij elke stap wordt een richting aangegeven voor het volgende aminozuur met waardes:
- 1: positieve **X**-richting
- -1: negatieve **X**-richting
- 2: positieve **Y**-richting
- -2: negatieve **Y**-richting
- 3: positieve **Z**-richting
- -3: negatieve **Z**-richting

### Scoresysteem
Hydrofobe aminozuren (H) hebben een voorkeur om naast elkaar te liggen, wat resulteert in stabiliserende "H-bonds." Polaire aminozuren (P) hebben geen voorkeur, terwijl Cysteïne (C) sterker bindt met andere C's en zwakker met H's. Het uiteindelijke doel is om eiwitten zo te vouwen dat de stabiliteitsscore zo laag mogelijk is. De score is op basis van het volgende:

- Een H-H binding: **-1**
- Een C-C binding: **-5**
- Een C-H binding: **-1**
- Geen binding: **0**

## Aanpak Algoritmen
In dit project zijn de volgende algoritmes geïmplementeerd:
1. Random: Genereert random vouwingen en selecteert het resultaat met de laagste score.
2. Greedy: Gebruikt een random algoritme maar met greedy keuzes om de vijf iteraties, waarbij het probeert de score te verlagen. Het selecteert uiteindelijk het resultaat met de laagste score.
3. Hill Climber: Genereert eerst een random vouwing. Vervolgens optimaliseert het de vouwing door random aanpassingen aan de vouwingsstructuur, waarbij vouwingen met een verlaging van de score of gelijk zijnde score worden geaccepteerd. Hiermee wordt de vouwing iteratief verbeterd.
4. Simulated Annealing: Gebruikt dezelfde methode als Hill Climber, alleen worden in het begin ook scores die hoger zijn dan de huidige score geaccepteerd met een afnemende waarschijnlijkheid die afhankelijk is van de temperatuur. Dit wordt gedaan om lokale minima te vermijden.

## Aan de slag
### Vereisten
Deze codebase is volledig geschreven in Python 3.12. In requirements.txt staan alle benodigde packages om de code succesvol te draaien. Deze zijn gemakkelijk te installeren via pip dmv. de volgende instructie:
```
pip install -r requirements.txt
```
Of via conda:
```
conda install --file requirements.txt
```

### Gebruik
Een voorbeeldje kan gerund worden door aanroepen van:
```python
python main.py
```
In de main.py script kunnen de argumenten van de `run`-functie worden aangepast om specifieke algoritmes of eiwitten te runnen. Als de argumenten `protein_sequence` en `algorithm` niet worden gegeven, worden alle eiwitten en algoritmes gerund. Als er een specifieke eiwit getest moet worden, kan er gebruik gemaakt worden van de variable `protein_sequences`. Dit bevat een lijst met alle mogelijke eiwitten voor deze case. Het selecteren van één van de eiwitten kan met de volgende voorbeelden:
```python
protein_sequence = protein_sequences[0]
protein_sequence = protein_sequences[-1]
```
Als er een specifieke algortime getst moet worden, kan er gebruik gemaakt worden van de parameter `algorithm` met de volgende voorbeelden:
```python
algorithm = "Greedy"
algorithm = "Hill Climber"
```
Of een plot wordt getoond en opgeslagen, kan worden gedaan met de parameters:
```python
show: bool = True
save: bool = True
```
In welke dimensie de vouwing plaatsvindt kan met behulp van de `dimension` parameter. Dit kan alleen maar `2` of `3` zijn:
```python
dimension = 2
dimension = 3
```
Het aantal herhalingen (repeats) en iteraties kan ook worden aangepast met behulp van:
```python
repeats = 10
iterations = 10000
```
Deze aanpassingen kunnen gecombineerd worden zoals:
```python
run(protein_sequence = protein_sequences[4], algorithm = "Simulated Annealing", show = True, save = False, dimension = 3, repeats = 1, iterations = 1000)
```
Een zeer uitgebreide documentatie kan gevonden worden in de docstring van de `run` functie met nog meer voorbeelden.

Om de algoritmes te vergelijken kan een boxplot worden geplot van de verschillende distributies van scores van alle algortimes. De `view`-functie kan de boxplots tonen en opslaan. Deze functie heeft ook een parameter om te kiezen welke eiwit getoond wordt voor de vergelijking. Als dit niet wordt meegegeven in de functie, zal het alle eiwitten tonen.
Ook in deze functie kan worden bepaald of het in 2D of 3D is en of de plot getoond en opgeslagen moet worden. Een voorbeeld van deze functie:
```python
view(protein_sequence = "all", dimension = 2, show_plot = True, save_plot = False)
```
Voor een uitgebreide documentatie van deze functie, kan de docstring bekeken worden.

### Brute Force
Om een optimale oplossing te vinden is er ook een brute force methode aangemaakt die alle mogelijk eiwitten bekijkt en beoordeeld. Dit is echter alleen toepasbaar op de eerste 2 eiwitten in 2D. Desalniettemin is er voor het idee een script geschreven om alle mogelijke vouwingen te genereren in `python` en in `C`, die in een `CSV` bestand worden opgeslagen. Vervolgens kan met een functie `brute_force` al deze vouwingen geëvalueerd worden. Om dit te runnen kan het codestukje van Brute Force uitgecommend worden om het te runnen. Om alle mogelijke vouwingen te genereren met de `C` implementatie, kan het volgende command gerund worden:
```
./generate_all_foldings HHPHHHPH
```
Hiervoor moet je wel in de `src/brute_force` directory zitten. De resultaten hiervan worden ook in de `brute_force` folder opgeslagen.
### Structuur
```python
Helix-Hunters/
├── main.py                # Hoofdscript
├── src/
│   ├── algorithms/        # Bevat alle algoritmen (Random, Greedy, etc.)
│   ├── brute_force/       # Bevat de brute force methodes
│   ├── classes/           # Bevat drie benodigde classes voor deze case
│   ├── utils/             # Hulpfuncties en constanten
│   ├── visualisation/     # Visualisatie functies
├── data/                  # Map voor gegenereerde data
├── results/               # Map voor gegenereerde resultaten
├── requirements.txt       # Vereistenbestand
├── README.md              # Documentatie
```

## Auteurs
- Leslie Wu
- Lara Ersoy