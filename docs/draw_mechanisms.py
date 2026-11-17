"""Render editable RM mechanism diagrams with Python 3, without dependencies."""
from html import escape
from pathlib import Path

OUT = Path(__file__).resolve().parent / "src/assets"


def svg(name, title, description, body, h=600):
    content = f'''<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="title desc" viewBox="0 0 900 {h}">
<title id="title">{escape(title)}</title><desc id="desc">{escape(description)}</desc>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#475569"/></marker></defs>
<style>text{{font-family:system-ui,-apple-system,sans-serif;fill:#172b4d;font-size:22px}}.small{{font-size:18px}}.head{{font-size:24px;font-weight:650}}.edge{{fill:none;stroke:#475569;stroke-width:2.2;marker-end:url(#arrow)}}.muted{{fill:#475569}}.mono{{font-family:ui-monospace,monospace;font-size:20px}}</style>
<rect width="900" height="{h}" rx="14" fill="#ffffff"/>
{body}
</svg>'''
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(content)


def text(x, y, s, cls="", anchor="middle"):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" class="{cls}">{escape(s)}</text>'


def box(x, y, w, h, title, lines=(), color="#f0f5fb"):
    body = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{color}" stroke="#94a3b8"/>'
    body += text(x + w / 2, y + 34, title, "head")
    for i, line in enumerate(lines):
        body += text(x + w / 2, y + 67 + i * 29, line)
    return body


def edge(path):
    return f'<path class="edge" d="{path}"/>'


def draw_model_state():
    body = [
        '<rect x="30" y="25" width="840" height="195" rx="12" fill="#f8fafc" stroke="#94a3b8"/>',
        text(450, 60, "ConstData · held by reference", "head"),
        box(55, 80, 370, 120, "ProblemDefinition",
            ("Fixed families and keys", "Relations and objective")),
        box(450, 80, 395, 120, "Study payload",
            ("Arrays, dicts, structs", "Read-only between rebuilds"), "#ecfdf5"),
        '<rect x="315" y="275" width="570" height="205" rx="12" fill="#f8fafc" stroke="#94a3b8"/>',
        text(650, 308, "LinearProblem", "head"),
        edge("M470 220V330"),
        text(490, 258, "Read for formulas", "small", "start"),
        box(30, 330, 245, 135, "Parameter point",
            ("Values may change", "Same concrete type", "isbits value"), "#f5f3ff"),
        edge("M275 397H335"),
        box(335, 330, 270, 135, "Private workspace",
            ("Bounds, matrix, RHS", "Objective + results")),
        edge("M605 397H655"),
        box(655, 330, 210, 135, "Backend",
            ("Persistent", "solver model"), "#ecfdf5"),
        box(30, 485, 400, 110, "Study payload changes",
            ("rebuild! before solving", "every sharing problem"), "#fffbeb"),
        box(455, 485, 415, 110, "Families or keys change",
            ("New definition + problem",), "#fffbeb"),
    ]
    return svg(
        "mechanism-model-state.svg",
        "Definition, study data and numerical state",
        "ConstData holds a fixed mathematical definition and study payload by "
        "reference. Keep the payload read-only between rebuilds. A changing "
        "parameter point is an isbits value with one fixed concrete type. The "
        "definition, payload and point feed a private numerical workspace and "
        "persistent solver backend owned by one LinearProblem. After changing "
        "study data, rebuild every sharing problem before solving. Changed "
        "families or keys require a new definition and problem.",
        "\n".join(body),
        h=600,
    )


def draw_solve_flow():
    """Show first construction, safe/selective evaluation, and outcome reuse."""
    body = text(450, 36, "From parameter point to solver outcome", "head")
    body += text(40, 75, "FIRST OUTCOME · construction and solve", "small muted", "start")
    body += box(40, 95, 250, 100, "Build numerical model", ["Assemble all values"])
    body += box(330, 95, 270, 100, "Initialize backend", ["Retain the model"])
    body += box(640, 95, 220, 100, "Optimize", ["Record outcome"], "#eaf5ed")
    body += edge("M290 145H330")
    body += edge("M600 145H640")
    body += '<path d="M40 223H860" stroke="#cbd5e1" stroke-width="1.5"/>'

    body += text(40, 255, "LATER CALLS · same concrete point type", "small muted", "start")
    body += box(40, 280, 370, 118, "Safe mode", [
        "Reevaluate all numerical formulas",
        "even at an identical point",
    ], "#fff4e6")
    body += box(490, 280, 370, 118, "Selective mode", [
        "Identical point: reuse outcome",
        "Distinct point: test conditions",
    ])
    body += box(490, 430, 370, 85, "Evaluate affected formulas", [
        "Conditions choose reevaluation",
    ])
    body += edge("M225 398V550")
    body += edge("M675 398V430")
    body += edge("M675 515V550")

    body += box(40, 550, 820, 118, "Compare numerical values with isequal", [
        "Send only actual differences through backend setters",
        "Retain objective-constant changes in public results",
    ], "#f0f0fb")
    body += edge("M450 668V688H225V710")
    body += edge("M450 688H675V710")
    body += box(40, 710, 370, 85, "Optimize", [
        "Numerical values changed",
    ], "#eaf5ed")
    body += box(490, 710, 370, 85, "Reuse outcome", [
        "Numerical values unchanged",
    ], "#eaf5ed")
    body += edge("M860 347H882V752H860")
    body += text(450, 829, "A reused nonoptimal outcome raises SolveError again.", "small muted")
    return svg(
        "mechanism-solve-flow.svg",
        "First construction and subsequent safe or selective solves",
        "Initial construction assembles numerical values and initializes the backend; the first "
        "solve optimizes. Rebuilding repeats construction before the next solve. Later "
        "safe-mode calls reevaluate every numerical "
        "formula, even at an identical point. Selective mode reuses the previous outcome at an "
        "identical point; at a distinct point, dependency conditions choose which formulas to "
        "reevaluate. Both update paths compare numerical values with isequal and send only "
        "actual native differences through backend setters. Objective-constant changes are "
        "retained even without a native coefficient setter. A changed numerical model is "
        "optimized; an unchanged model reuses its previous outcome. Reusing a nonoptimal "
        "outcome raises SolveError again. Retaining the backend model does not guarantee basis reuse.",
        body,
        h=850,
    )


def draw_result_lifetime():
    body = text(450, 36, "One live solution, retained snapshots", "head")
    events = [(40, "solve(rhs=2)", "optimal", "#f0f5fb"),
              (250, "solve(rhs=4)", "optimal", "#f0f5fb"),
              (460, "solve(rhs=20)", "SolveError", "#fff1f2"),
              (670, "solve(rhs=3)", "optimal", "#f0f5fb")]
    for x, title, outcome, color in events:
        body += box(x, 70, 190, 85, title, (outcome,), color)
    body += text(40, 191, "Live Solution · the same object", "head", "start")
    for x, value, detail, color in [(40, "2.0", "first result", "#ecfdf5"),
                                    (250, "4.0", "refreshed", "#ecfdf5"),
                                    (460, "invalid", "queries raise", "#fff1f2"),
                                    (670, "3.0", "refreshed", "#ecfdf5")]:
        body += box(x, 210, 190, 95, value, (detail,), color)
    for x in (230, 440, 650):
        body += edge(f"M{x} 257H{x+20}")
    body += text(40, 350, "Values retained from the first solve", "head", "start")
    for x in (40, 250, 460, 670):
        body += box(x, 375, 190, 95, "2.0", ("scalar + snapshot",), "#f5f3ff")
    body += text(450, 519, "saved_value and saved_variables keep 2.0.")
    body += text(450, 554, "After recovery, recovered === solution.")
    svg("mechanism-result-lifetime.svg", "Live solution and retained snapshots",
        "The executed example minimizes x with zero to ten bounds. At rhs=2 the live "
        "solution, saved scalar and materialized variable snapshot have value two. "
        "At rhs=4 the same live solution reports four; the retained scalar and snapshot "
        "stay two. At rhs=20 the solve is infeasible, so live queries raise while the "
        "retained values stay two. Recovery at rhs=3 refreshes the same solution object "
        "to three; retained values are still two.", body, h=585)


if __name__ == "__main__":
    draw_model_state()
    draw_solve_flow()
    draw_result_lifetime()
