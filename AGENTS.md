# AGENTS.md

Instructions for coding agents that work in this repository.

## Coding conventions

### Line endings

- Use Unix-style line endings (LF, `\n`) in all text files.
- Do not use Windows-style line endings (CRLF, `\r\n`).
- End each text file with one newline character.
- Some legacy files still use CRLF (for example `hemoFlow.cpp`, `opening.cpp`,
  `template_config.xml`). To list them, run `git ls-files --eol | grep crlf`.
- To convert a file to LF, run `sed -i 's/\r$//' <file>`.
- Do a line-ending conversion in a separate commit.
  Do not mix it with code changes, because it makes the diff hard to read.
