# Scientific Data Integrity Policy

**CRITICAL RULE: DO NOT FALSIFY OR MOCK SCIENTIFIC DATA**

This project enforces a strict zero-tolerance policy for data falsification, mocking, or using hardcoded placeholder values in analytical pipelines.

As an AI agent working on this repository, you must adhere to the following rules:

1. **No Hardcoded Metrics:** Never hardcode values that should be computed from data (e.g., proximity scores, p-values, fold changes, cell counts). If a module is incomplete, use `raise NotImplementedError()` to halt execution.
2. **No Random Data Fallbacks:** Never use `numpy.random` (or similar tools) to generate "highly realistic mock data" when an upstream step fails or data is missing. If data cannot be found or parsed, the pipeline MUST crash with a clear, descriptive Error.
3. **No Placeholders in Output DataFrames:** If building a pipeline structure before the computation logic is ready, do not populate output CSVs or DataFrames with fake rows just to "make it work." Leave the DataFrame empty or halt execution.
4. **Data Provenance:** All gene sets, signatures, and constants must be traceable to their source (e.g., centralized config files like `druggability_config.py` or cited literature). Do not invent biological mechanisms.

**Failure to follow these rules damages the credibility of the human scientist and the integrity of the research.** If you encounter existing code that violates these rules, your immediate priority is to rip it out and replace it with hard errors, then notify the user.
