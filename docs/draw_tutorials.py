"""Draw editable SVG tutorial schematics; run with Python 3, no dependencies."""
from pathlib import Path
from html import escape
out=Path(__file__).resolve().parent / 'src/assets'

def svg(name,title,description,body,h=360):
 s=f'''<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="title desc" viewBox="0 0 1000 {h}">
<title id="title">{escape(title)}</title><desc id="desc">{escape(description)}</desc>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#475569"/></marker></defs>
<style>text{{font-family:system-ui,-apple-system,sans-serif;fill:#172b4d;font-size:18px}}.small{{font-size:15px}}.head{{font-size:21px;font-weight:600}}.edge{{fill:none;stroke:#475569;stroke-width:2;marker-end:url(#arrow)}}</style>
<rect width="1000" height="{h}" rx="12" fill="#ffffff"/>
{body}</svg>'''
 (out/name).write_text(s)

def text(x,y,s,cls='',anchor='middle'):
 return f'<text x="{x}" y="{y}" text-anchor="{anchor}" class="{cls}">{escape(s)}</text>'
def box(x,y,w,h,title,lines=(),color='#f0f5fb'):
 b=f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" fill="{color}" stroke="#94a3b8"/>'
 b+=text(x+w/2,y+30,title,'head')
 for i,line in enumerate(lines): b+=text(x+w/2,y+60+i*25,line)
 return b
def edge(d):return f'<path class="edge" d="{d}"/>'

b=text(500,32,'Two-generator dispatch: one-hour periods','head')
b+=box(35,65,280,100,'Generator 1',('Capacity 80 MW; cost 40/MWh',))
b+=box(35,205,280,100,'Generator 2',('Capacity 60 MW; cost 70/MWh',))
b+=box(410,130,180,100,'Balance',('p₁ + p₂ = Dₜ',))
b+=box(710,100,255,155,'Demand',('t = 1: 100 MW','t = 2: 110 MW'))
b+=edge('M315 115 L375 115 L410 160')+text(375,100,'p₁','small')
b+=edge('M315 255 L375 255 L410 200')+text(370,280,'p₂','small')
b+=edge('M590 180 L710 180')
b+=text(650,155,'Supply','small')
b+=text(500,340,'Solution: (80, 20) MW, then (80, 30) MW; costs 4600, then 5300.')
svg('tutorial-dispatch.svg','Dispatch model','Two bounded generators supply a changing demand. Generator 1 costs 40 per MWh and has capacity 80 MW. Generator 2 costs 70 and has capacity 60. Demand rises from 100 to 110 MW.',b,365)

b=text(500,30,'Network flow at t = 1','head')
b+=edge('M220 200 L450 140')+text(310,115,'Flow 10; capacity 10','small')+text(310,139,'Cost 2/unit','small')
b+=edge('M550 140 L780 200')+text(690,115,'Flow 7; capacity 8','small')+text(690,139,'Cost 1/unit','small')
b+=edge('M165 281 L165 335 L835 335 L835 281')+text(500,320,'Direct arc: flow 0; capacity 4; cost 4/unit','small')
for x,y,title,sub in [(165,220,'North','Injection +10'),(500,125,'Centre','Injection −3'),(835,220,'South','Injection −7')]:
 b+=f'<circle cx="{x}" cy="{y}" r="61" fill="#f0f5fb" stroke="#94a3b8"/>'+text(x,y-5,title,'head')+text(x,y+23,sub,'small')
b+=text(500,390,'Conservation: outgoing − incoming = injection. Total cost: 10 × 2 + 7 × 1 = 27.')
svg('tutorial-network.svg','Minimum-cost network flow','At the first time step North supplies 10 units, Centre consumes 3 and South consumes 7. Ten units travel from North to Centre and seven from Centre to South. No flow uses the direct arc.',b,415)

b=text(500,30,'A reservoir transition at period t','head')
b+=box(340,105,280,130,'Reservoir',('0 ≤ sₜ ≤ capacity','Water balance'))
b+=edge('M35 170 L340 170')+text(180,145,'Incoming stock sₜ₋₁')
b+=edge('M480 63 L480 105')+text(480,56,'Inflow Δₜ aωₜ')
b+=edge('M620 140 L760 140')+box(765,90,210,90,'Turbine',('Revenue Δₜ πωₜ ηₜ uₜ',))
b+=edge('M620 208 L760 233')+text(825,235,'Spill Δₜ zₜ')
b+=edge('M480 235 L480 305')+text(345,278,'Outgoing stock sₜ')
b+=box(245,310,485,68,'Next period, or terminal value',(), '#f6f4fb')
b+=text(500,400,'sₜ − sₜ₋₁ + Δₜ(uₜ + zₜ) = Δₜ aωₜ. At T: α ≤ βₖ sT + γₖ for every cut k.','small')
svg('tutorial-reservoir.svg','Reservoir stage balance','Incoming stock and inflow are split between outgoing stock, turbine release and spill. At the horizon the remaining stock enters supplied terminal-value cuts.',b,425)

b=text(500,30,'SDDP: sample states, then improve future-cost cuts','head')
for x,t in [(40,1),(375,2),(710,3)]:
 b+=box(x,90,250,105,f'Stage {t}',('Solve LP at observed inflow', 'No cost after stage 3' if t == 3 else 'Pass outgoing stock forward'))
b+=edge('M290 140 L375 140')+edge('M625 140 L710 140')
b+=text(500,66,'Forward: sample one inflow path; no knowledge of future outcomes.','small')
b+=edge('M840 195 L840 280 L625 280')+edge('M625 280 L290 280')+edge('M290 280 L165 280 L165 195')
b+=text(500,265,'Backward: enumerate both inflows, average values and duals.','small')
b+=text(500,320,'Insert cuts in the preceding LP, rebuild after cut-data mutation, repeat.','small')
b+=text(500,355,'Check: lower bound ≤ exact scenario-tree optimum ≤ expected feasible-policy cost.','small')
svg('tutorial-sddp.svg','Forward and backward SDDP passes','Three stage LPs pass stock forward along a sampled inflow path. The backward pass enumerates the two inflows at each visited state, averages values and balance duals, then updates future-cost cuts in the preceding stage.',b,380)
