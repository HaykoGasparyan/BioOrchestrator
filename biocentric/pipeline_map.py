"""Compact horizontal pipeline map, driven exclusively by recorded evaluation state."""
import html
import json
import textwrap

COLORS = {'gray': '#c4ceda', 'yellow': '#e4b34c', 'green': '#32aa87'}
STAGES = [('structure_preparation', 'Protein', '#7781d9'), ('pocket_identification', 'Pocket', '#34a8b0'),
          ('docking', 'Docking', '#6689d5'), ('scoring', 'Scoring', '#a47ac7'), ('ligands', 'Ligands', '#5cae9c')]


def graph_data(selected, labels, review=None, state=None, comparisons=(), locked=None, identity=None):
    state, review = state or {}, review or {}
    identity = review.get('identity') or identity or {}
    target = (identity.get('gene_names') or [identity.get('accession') or 'Protein'])[0]
    nodes, stages = [], []
    for phase, title, accent in STAGES:
        if phase == 'structure_preparation':
            group = [{'id': 'protein', 'tool': 'protein', 'label': target,
                      'detail': review.get('protocol', {}).get('preparation', 'Protein structure from the reviewed experimental reference')}]
        elif phase == 'pocket_identification':
            if review.get('candidates') and not any('pocket' in p for p in review['candidates']):
                group = [{'id': 'pocket', 'tool': 'pocket', 'label': 'Crystal pocket',
                          'detail': review.get('protocol', {}).get('pocket', 'Legacy crystal-defined site')}]
            else:
                group = [{'id': phase + ':' + item['tool'], 'tool': item['tool'],
                          'label': labels.get(item['tool'], item['tool']),
                          'detail': json.dumps(item.get('parameters', {}), ensure_ascii=False)}
                         for item in selected.get(phase, [])]
        elif phase == 'ligands':
            group = [{'id': 'ligands', 'tool': 'ligands', 'label': 'Evaluation ligands', 'detail': 'Measured compounds and the crystal control in the frozen evaluation set'}]
        else:
            group = [{'id': phase + ':' + item['tool'], 'tool': item['tool'], 'label': labels.get(item['tool'], item['tool']),
                      'detail': json.dumps(item.get('parameters', {}), ensure_ascii=False)} for item in selected.get(phase, [])]
        for node in group:
            node.update(phase=phase, accent=accent)
        stages.append(group)
        nodes.extend(group)
    paths = {p['id']: ['protein', 'pocket_identification:' + p['pocket']['tool'] if 'pocket' in p else 'pocket', 'docking:' + p['docking']['tool'], 'scoring:' + p['scoring']['tool'], 'ligands'] for p in review.get('candidates', [])}
    scores = {r['id']: r for r in comparisons}
    active = paths.get(state.get('active_combination'), []) if state.get('status') == 'running' else []
    final = paths.get((locked or {}).get('candidate', {}).get('id'), [])
    edges = []
    for left, right in zip(stages, stages[1:]):
        for a in left:
            for b in right:
                covered = [key for key, path in paths.items() if any(x == a['id'] and y == b['id'] for x, y in zip(path, path[1:]))]
                values = [scores[key]['color'] for key in covered if key in scores]
                color = 'green' if values and len(values) == len(covered) and all(v == 'green' for v in values) else 'yellow' if any(v != 'gray' for v in values) else 'gray'
                detail = f"{a['label']} → {b['label']} · {color}"
                if covered:
                    detail += f' · {len(covered)} candidate path(s)'
                    for key in covered:
                        r = scores.get(key)
                        if r:
                            detail += f" | {r.get('pocket', 'reference_box')} → {r['docking']} → {r['scoring']}: {r['completed']}/{r['expected']} compounds; ρ=" + (f"{r['rho']:.3f}" if r['rho'] is not None else 'unavailable')
                            if r['rho'] is not None:
                                detail += f" (correlation n={r.get('ranked_measurements', '?')}/{r['expected']})"
                else:
                    detail += ' · Not evaluated'
                edges.append({'source': a['id'], 'target': b['id'], 'color': color,
                              'running': any(x == a['id'] and y == b['id'] for x, y in zip(active, active[1:])),
                              'fixed': any(x == a['id'] and y == b['id'] for x, y in zip(final, final[1:])), 'detail': detail})
    current_phase = {'preparation': 'structure_preparation', 'pocket': 'pocket_identification', 'validation': 'ligands'}.get(state.get('stage'), state.get('stage'))
    for n in nodes:
        n.update(running=n['id'] in active, current=n['id'] in active and n['phase'] == current_phase, fixed=n['id'] in final)
    compounds = ([{'compound_id': 'crystal:' + review['reference']['component'], 'endpoint': None}] if review.get('reference', {}).get('component') else []) + review.get('records', [])
    ligand_rows = []
    seen = set()
    for record in compounds:
        if record['compound_id'] in seen:
            continue
        seen.add(record['compound_id'])
        active_ligand = state.get('status') == 'running' and state.get('active_compound') == record['compound_id']
        ligand_rows.append({'id': record['compound_id'], 'active': active_ligand,
                            'detail': record['compound_id']})
    return {'nodes': nodes, 'stages': stages, 'edges': edges, 'ligands': ligand_rows}


def graph_height(data):
    return max(290, max(len(data['stages'][1]), len(data['stages'][2]), len(data['stages'][3]), 1) * 62 + 100)


def graph_html(data):
    escape = lambda value: html.escape(str(value), quote=True)
    height = graph_height(data)
    positions, fragments = {}, []
    centers = [68, 210, 360, 515, 674]
    middle = height / 2 + 8
    for index, (group, (_, title, accent)) in enumerate(zip(data['stages'], STAGES)):
        x = centers[index]
        fragments.append(f'<text x="{x}" y="28" text-anchor="middle" class="stage" fill="{accent}">{escape(title)}</text>')
        fragments.append(f'<line x1="{x-20}" x2="{x+20}" y1="39" y2="39" stroke="{accent}" opacity=".2"/>')
        for i, node in enumerate(group):
            positions[node['id']] = (x, middle + (i - (len(group)-1)/2) * 62)
    for edge in data['edges']:
        x,y = positions[edge['source']]; xx,yy = positions[edge['target']]
        start,end = x+49,xx-49
        path = f'M{start},{y} C{start+26},{y} {end-26},{yy} {end},{yy}'
        attrs = f'data-graph-key="{escape("edge:" + edge["source"] + "|" + edge["target"])}" tabindex="0" role="button" aria-label="{escape(edge["detail"])}" data-detail="{escape(edge["detail"])}"'
        fragments.append(f'<g {attrs}><path d="{path}" fill="none" stroke="transparent" stroke-width="14"/>')
        edge_color = '#15966d' if edge['fixed'] else COLORS[edge['color']]
        edge_class = 'edge selected-edge' if edge['fixed'] else 'edge'
        if edge['fixed']:
            fragments.append(f'<path d="{path}" class="fixed-edge"/>')
        fragments.append(f'<path d="{path}" class="{edge_class}" stroke="{edge_color}"/>')
        # A visible connection dot also gives straight SVG paths a nonzero hit box.
        fragments.append(f'<circle cx="{(start+end)/2}" cy="{(y+yy)/2}" r="2.5" fill="{edge_color}"/>')
        fragments.append(f'<path d="{path}" class="flow" style="opacity:{1 if edge["running"] else 0}"/>')
        fragments.append(f'<title>{escape(edge["detail"])}</title></g>')
    for node in data['nodes']:
        x,y = positions[node['id']]
        detail = node['label'] + ' · ' + node['detail']
        classes = 'node' + (' current' if node['current'] else '') + (' fixed' if node['fixed'] else '')
        fragments.append(f'<g data-graph-key="{escape("node:" + node["id"])}" class="{classes}" tabindex="0" role="button" aria-label="{escape(detail)}" data-detail="{escape(detail)}" style="--accent:{node["accent"]}">')
        if node['id'] in {'protein', 'pocket'}:
            fragments.append(f'<rect x="{x-46}" y="{y-51}" width="92" height="102" rx="21"/>')
            if node['id'] == 'protein':
                fragments.append(f'<g transform="translate({x},{y-15})" fill="none" stroke-linecap="round"><path d="M-22,-18 C24,-31 -26,30 22,19 M-15,-24 C-35,0 35,3 15,27" stroke="#8c90df" stroke-width="5"/><path d="M-19,-9 L12,-16 M-10,1 L9,-3 M-10,15 L17,7" stroke="#70c4cb" stroke-width="3"/></g>')
            else:
                fragments.append(f'<g transform="translate({x},{y-16})"><path d="M-27,-7 Q-10,-33 11,-21 L24,-4 L15,20 L-12,22 L-26,7 Z" fill="#e2f6f2" stroke="#64bfc0" stroke-width="2"/><circle r="11" fill="#fff" stroke="#54b5ad" stroke-width="1.5" stroke-dasharray="3 3"/><circle r="4" fill="#63b9b4"/></g>')
            fragments.append(f'<text x="{x}" y="{y+30}" text-anchor="middle" class="name">{escape(node["label"][:18])}</text>')
            if node.get('subtitle'):
                fragments.append(f'<text x="{x}" y="{y+44}" text-anchor="middle" class="sub">{escape(node["subtitle"])}</text>')
        elif node['id'] == 'ligands':
            ligands = data['ligands']
            visible = ligands[:7]
            active_ligand = next((r for r in ligands[7:] if r['active']), None)
            if active_ligand:
                visible[-1] = active_ligand
            box_height = max(64, len(visible)*20 + 29)
            top = y - box_height/2
            fragments.append(f'<rect x="{x-49}" y="{top}" width="98" height="{box_height}" rx="15"/>')
            if not ligands:
                fragments.append(f'<text x="{x}" y="{y+4}" text-anchor="middle" class="sub">Awaiting test set</text>')
            for i,ligand in enumerate(visible):
                ly = top+21+i*20
                fragments.append(f'<g data-graph-key="{escape("ligand:" + ligand["id"])}" tabindex="0" role="button" aria-label="{escape(ligand["detail"])}" data-detail="{escape(ligand["detail"])}"><circle cx="{x-37}" cy="{ly-3}" r="2.5" fill="{"#04b8c4" if ligand["active"] else "#b4c2cb"}"/><text x="{x-28}" y="{ly}" class="ligand">{escape(ligand["id"])}</text><title>{escape(ligand["detail"])}</title></g>')
            if len(ligands)>7:
                shown = {ligand['id'] for ligand in visible}
                remaining = ' · '.join(ligand['id'] for ligand in ligands if ligand['id'] not in shown)
                fragments.append(f'<text x="{x}" y="{top+box_height-9}" text-anchor="middle" class="sub"><title>{escape(remaining)}</title>+{len(ligands)-7} ligands</text>')
        else:
            fragments.append(f'<rect x="{x-49}" y="{y-23}" width="98" height="46" rx="13"/>')
            lines = textwrap.wrap(node['label'], 15)[:2]
            for j,line in enumerate(lines):
                fragments.append(f'<text x="{x}" y="{y+4+(j-(len(lines)-1)/2)*14}" text-anchor="middle" class="name">{escape(line)}</text>')
        fragments.append(f'<title>{escape(detail)}</title></g>')
    return '''<!doctype html><html><head><meta charset="utf-8"><style>
*{box-sizing:border-box}body{margin:0;color:#3b506c;font:12px system-ui;background:transparent}.board{border:1px solid #e4eaf1;border-radius:22px;background:radial-gradient(ellipse at 5% 50%,#e7ecff99,transparent 45%),radial-gradient(ellipse at 95% 40%,#ddf4ef88,transparent 50%),#fcfdff;overflow:auto;padding:8px 0}svg{width:100%;min-width:420px;display:block}.stage{font-size:13px;font-weight:650;letter-spacing:.2px}.edge{transition:stroke .25s ease;fill:none;stroke-width:1.6;stroke-dasharray:2 5;stroke-linecap:round}.fixed-edge{fill:none;stroke:#15966d;stroke-width:7;opacity:.2}.selected-edge{stroke-width:2.8;stroke-dasharray:none}.flow{fill:none;stroke:#04b8c4;stroke-width:2.4;transition:opacity .25s ease;filter:drop-shadow(0 0 3px #52d6df)}.node>rect{transition:stroke .25s ease,stroke-opacity .25s ease,fill .25s ease,filter .25s ease;fill:#fff;stroke:var(--accent);stroke-opacity:.24;stroke-width:1.3;filter:drop-shadow(0 4px 5px #30477708)}.node.current>rect{stroke:#04b8c4;stroke-opacity:1;stroke-width:2;filter:drop-shadow(0 0 5px #54d6df66)}.node.fixed>rect{fill:#eefbf4;stroke:#15966d;stroke-opacity:1;stroke-width:2}.name{font-size:12px;font-weight:550;fill:#40516a}.sub{font-size:9px;fill:#8997aa}.ligand{font-size:9px;fill:#718698}g[role=button]{cursor:pointer;outline:none}g[role=button]:focus>.edge,g[role=button]:hover>.edge{stroke-width:3}g[role=button]:focus>rect,g[role=button]:hover>rect{stroke-opacity:.85;stroke-width:2}.legend{display:flex;gap:12px;flex-wrap:wrap;padding:15px 5px 8px;font-size:10px;color:#8796a7}.dot{display:inline-block;width:6px;height:6px;border-radius:50%;margin-right:4px}#detail{margin:0 5px;line-height:1.6;color:#77899f;font-size:10px;max-height:62px;overflow:auto}@media(prefers-reduced-motion:reduce){.flow,.edge,.node>rect{transition:none}}
</style></head><body><div class="board"><svg role="group" aria-label="Horizontal pipeline connections" viewBox="0 0 742 ''' + str(height) + '''">''' + ''.join(fragments) + '''</svg></div><div class="legend"><span><i class="dot" style="background:#c4ceda"></i>Untested</span><span><i class="dot" style="background:#e4b34c"></i>Limited</span><span><i class="dot" style="background:#32aa87"></i>Supported</span><span><i class="dot" style="background:#04b8c4"></i>Running</span><span><i class="dot" style="background:#15966d"></i>Locked</span></div><p id="detail"></p><script>document.querySelectorAll('[data-detail]').forEach(el=>{const show=e=>{e.stopPropagation();document.getElementById('detail').textContent=el.dataset.detail};el.addEventListener('click',show);el.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();show(e)}})});</script></body></html>'''
