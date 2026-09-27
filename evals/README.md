# Model evaluation boundary

The application has a separate optional model adapter (`biocentric/agent.py`). It does not share scientific grading with the language model.

`tests/test_agent.py` is an **offline software-contract suite**, using an injected fake Responses client. It checks schema constraints, unknown/missing citations, fabricated computational-result categories, unknown tool names, untrusted evidence placement, local-only privacy, structured response parsing, call traces and unread-evidence citation rejection. Passing these tests does **not** establish the scientific accuracy of a live model.

The live evaluation harness can be run with a configured API key and an explicitly chosen model:

```sh
.venv/bin/python -m biocentric.evaluate_agent --model YOUR_MODEL --output runs/model-evaluation.json
```

It runs bounded public, synthetic evaluation cases for endpoint fidelity, absence of completed calculations, target ambiguity and retrieved prompt injection. The report saves actual structured outputs, trace events and deterministic rubric flags for scientist review. Flags detect obvious failures; they are not full entailment or causal-reasoning graders. Review source support manually before using a model for research recommendations. Token use is bounded by 3 calls × 900 output tokens per case; API charges apply.

No live model evaluation was performed during initial delivery because `OPENAI_API_KEY` was not configured. No model scientific-performance metric is claimed. Scientific tool behavior is tested separately, and real local docking results are unrelated to LLM benchmark performance.
