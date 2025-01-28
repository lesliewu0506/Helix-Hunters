# Protein Pow(d)er
## Introductie case
Eiwitten zijn essentiële moleculen in het menselijk lichaam, opgebouwd uit lange ketens van aminozuren. De specifieke vouwing van deze eiwitten bepaalt hun functie. Verkeerd gevouwen eiwitten kunnen leiden tot ernstige ziekten zoals kanker, Alzheimer en taaislijmziekte. Het doel van dit project is om algoritmes te maken die de meest stabiele vouwing vinden voor eitwitstructuren.

### Vouwing
Dit project maakt gebruik van een 3D-grid (optioneel ook 2D) om eiwitten te modelleren.Hierbij wordt elke aminozuur stap voor stap geplaatst in de grid. Bij elke stap wordt een richting aangegeven voor het volgende aminozuur met waardes:
- 1: positieve **X**-richting
- -1: negatieve **X**-richting
- 2: positieve **Y**-richting
- -2: negatieve **Y**-richting
- 3: positieve **Z**-richting
- -3: negatieve **Z**-richting
- 0: einde van eiwitketen, geen buigingen meer

Elke vouwing van een aminozuur heeft een hoek van 90 graden. Een lijst houdt deze waardes bij, een voorbeeld hiervan ziet er zo uit:
```python
amino_directions: = [-2, 1, -3, 2, 3, 1, -2, 3, -2, -1, 2, -1, 2, 1, 3, 1, -3, 2, -3, 0]
```
### Scoresysteem
Hydrofobe aminozuren (H) hebben een voorkeur om naast elkaar te liggen, dat resulteert in stabiliserende "H-bond". Polaire aminozuren (P) hebben geen voorkeur, terwijl Cysteïne (C) sterker bindt met andere C's en zwakker met H's. Het uiteindelijke doel is om eiwitten zo te vouwen dat de stabiliteitsscore zo laag mogelijk is. De score is op basis van het volgende:

- Een H-H binding: **-1**
- Een C-C binding: **-5**
- Een C-H binding: **-1**
- Overige binding: **0**

## Aanpak Algoritmes
In dit project zijn de volgende algoritmes geïmplementeerd:
1. Random: Genereert random vouwingen en slaat het resultaat met de laagste score op.
2. Greedy: Gebruikt een random algoritme maar met greedy keuzes om de vijf iteraties, waarbij het probeert de score te verlagen. Het selecteert uiteindelijk het resultaat met de laagste score en slaat het op.
3. Hill Climber: Genereert eerst een random vouwing. Vervolgens optimaliseert het de vouwing door random aanpassingen aan de vouwingsstructuur, waarbij vouwingen met een verlaging van de score of gelijk zijnde scores worden geaccepteerd. Hiermee wordt de vouwing iteratief verbeterd. De vouwing met de laagste score wordt opgeslagen.
4. Simulated Annealing: Gebruikt dezelfde methode als Hill Climber, alleen worden in het begin ook scores die hoger zijn dan de huidige score geaccepteerd met een afnemende waarschijnlijkheid die afhankelijk is van de temperatuur. Dit wordt gedaan om lokale minima te vermijden. Met trail en error zijn de volgende waardes bepaalt voor Simulated Annealing: 
    - Temperatuur: 4
    - Exponentiële temperatuur functie: $T_{current} = T_{start} * (0.999)^{iterations}$
    - Acceptatie functie: $p = e^{-(rating_{new} - rating_{old}) / T}$

## Aan de slag
### Vereisten
Deze codebase is volledig geschreven in Python 3.12. In `requirements.txt` staan alle benodigde packages om de code succesvol te draaien. Deze zijn gemakkelijk te installeren via pip dmv. de volgende instructie:
```
pip install -r requirements.txt
```
Of via conda:
```
conda install --file requirements.txt
```
### Ubuntu WSL Problemen
Het kan voorkomen dat in een ubuntu terminal de volgende message wordt weergegeven:
```
A module that was compiled using NumPy 1.x cannot be run in
NumPy 2.2.1 as it may crash. To support both 1.x and 2.x
versions of NumPy, modules must be compiled with NumPy 2.0.
Some module may need to rebuild instead e.g. with 'pybind11>=2.12'.

If you are a user of the module, the easiest solution will be to
downgrade to 'numpy<2' or try to upgrade the affected module.
We expect that some modules will need time to support NumPy 2.
```
Een oplossing hiervoor is om een virtual environment te maken en daarin een oudere versie van numpy te installeren. Het runnen van de commands in command prompt of powershell zou geen problemen moeten veroorzaken.

### Run Functie
Met de CLI kan de `run` command worden gebruikt om de algoritmes te runnen voor de eiwitten. Als de geen opties worden meegegeven, dan worden alle algoritmes gerund voor alle eiwitten met 1 repeat 10000 iteraties. Een voorbeeld ziet er zo uit:
```
$ python main.py run
```
- Als er een specifieke eiwit getest moet worden, kan er gebruik gemaakt worden van de optie `-p/--protein`:
```
$ python main.py run -p HHPHHHPHPHHHPH
$ python main.py run --protein HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH
```
- Als er een specifieke algortime getst moet worden, kan er gebruik gemaakt worden van de optie `-a/--algorithm`:
```
$ python main.py run -a Random
$ python main.py run -a Greedy
$ python main.py run --algorithm "Hill Climber"
$ python main.py run --algorithm "Simulated Annealing"
```
- Of een plot wordt getoond en opgeslagen, kan worden gedaan met de volgende flags:
```python
$ python main.py run --graph # Plotten
$ python main.py run --save  # Opslaan
```
- In welke dimensie de vouwing plaatsvindt, kan met behulp van de `-d/--dimension` optie:
```
$ python main.py run -d 3
$ python main.py run --dimension 2
```
- Het aantal herhalingen (repeats) `-r/--repeats` en iteraties `-i/--iterations` kan ook worden aangepast met behulp van:
```
$ python main.py run -r 5
$ python main.py run -iterations 1000
```
Deze aanpassingen kunnen gecombineerd worden zoals:
```
$ python main.py run -p HHPHHHPHPHHHPH -a "Hill Climber" -d 2 --graph -r 5 -i 1000 --graph --save
```
Een zeer uitgebreide documentatie kan gevonden worden in de docstring van de `run` functie met nog meer voorbeelden.

### View Functie

Om de algoritmes te vergelijken kan een boxplot worden geplot van de verschillende distributies van scores van alle algortimes. De `view` command kan de boxplots tonen en opslaan. Deze command heeft ook een optie om te kiezen welke eiwit getoond wordt voor de vergelijking. Als dit niet wordt meegegeven in de functie, zal het alle eiwitten tonen.
Ook in deze functie kan worden bepaald of het in 2D of 3D is en of de plot getoond en opgeslagen moet worden met dezelfde opties als de `run` command:
```
$ python main.py view --save
$ python main.py view -p HPHPPHHPHPPHPHHPPHPH -d 2 --graph --save
```
Voor een uitgebreide documentatie van deze functie, kan de docstring bekeken worden.

### Brute Force
Om een optimale oplossing te vinden is er ook een brute force methode aangemaakt die alle mogelijk eiwitten bekijkt en beoordeeld. Dit is echter alleen toepasbaar op de eerste 2 eiwitten in 2D. Desalniettemin is er voor het idee een script geschreven om alle mogelijke vouwingen te genereren in `python` en in `C`, die in een `CSV` bestand worden opgeslagen. Vervolgens kan met een functie `brute_force` al deze vouwingen geëvalueerd worden. Om deze functies te runnen, kan het command `bruteforce` geroepen worden:
```
$ python main.py bruteforce
```

Om alle mogelijke vouwingen te genereren met de `C` implementatie, kan het volgende command gerund worden:
```
./generate_all_foldings HHPHHHPH
```
Hiervoor moet je wel in de `src/brute_force` directory zitten. De resultaten hiervan worden ook in de `brute_force` folder opgeslagen.

### Structuur
```python
Helix-Hunters/
├── main.py                # Hoofdscript
├── src/
│   ├── algorithms/        # Bevat alle algoritmes (Random, Greedy, etc.)
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