## Lumrik should favor **ownership-driven, concise Rust APIs**.

Put data and the behaviour that naturally belongs to that data together. A type should normally implement operations on the state it owns instead of requiring unrelated callers to inspect its internals, reconstruct its logic, or coordinate several low-level functions.

Keep public APIs small and intention-oriented. Hide implementation machinery inside the crate that owns it. Callers should say **what they want done**, not reproduce **how it is done**.

Do not apply this dogmatically. Avoid Java-style getter/setter pairs, excessive wrapper types, unnecessary traits, and abstraction for abstraction's sake. Direct field access and free functions are perfectly appropriate when they make the code simpler and ownership remains clear.

Prefer existing objects as sources of truth. When information already exists in an object, provide concise operations that derive or consume it rather than copying its fields into parallel state.

Avoid duplicated knowledge: filenames, lifecycle rules, serialization conventions, state transitions, and other invariants should have one owner and one implementation.

A good Lumrik API should make the normal call site **short, obvious, and difficult to misuse**.

When refactoring, inspect the actual implementation and its callers first. Preserve computational behaviour unless intentionally changing it, but do not preserve accidental APIs or AI-generated architectural baggage merely for backward compatibility.

The goal is pragmatic, readable Rust:

**own the data → own its natural behaviour → expose the smallest useful interface.**

