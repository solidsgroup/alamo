from pygls.server import LanguageServer
from lsprotocol.types import (
    CompletionItem, CompletionList, CompletionParams, Hover, HoverParams, MarkupContent, MarkupKind, Diagnostic, DiagnosticSeverity, Position, Range, TextDocumentPositionParams
)
import json

# Load input key definitions from a JSON file
with open("input_schema.json", "r") as f:
    INPUT_SCHEMA = json.load(f)

server = LanguageServer("InputLSP", "0.1")

def get_key_info(key: str):
    """Retrieve information about a key from the schema."""
    return INPUT_SCHEMA.get(key, None)

@server.feature("textDocument/completion")
def completions(params: CompletionParams):
    """Provide autocompletion for input keys."""
    return CompletionList(is_incomplete=False, items=[
        CompletionItem(label=key, detail=value.get("description", ""))
        for key, value in INPUT_SCHEMA.items()
    ])

@server.feature("textDocument/hover")
def hover(params: HoverParams):
    """Provide hover documentation for known keys."""
    position = params.position
    key = ""  # Logic to extract key from position (simplified for now)
    key_info = get_key_info(key)
    
    if key_info:
        return Hover(contents=MarkupContent(kind=MarkupKind.Markdown, value=f"**{key}**\n\n{key_info['description']}"))
    return None

@server.feature("textDocument/publishDiagnostics")
def validate_document(params):
    """Validate input file and highlight unknown keys."""
    uri = params.text_document.uri
    text = server.workspace.get_document(uri).source
    diagnostics = []
    
    for i, line in enumerate(text.split("\n")):
        key = line.split("=")[0].strip()
        if key and key not in INPUT_SCHEMA:
            diagnostics.append(Diagnostic(
                range=Range(start=Position(line=i, character=0), end=Position(line=i, character=len(key))),
                message=f"Unknown input key: {key}",
                severity=DiagnosticSeverity.Warning
            ))
    
    server.publish_diagnostics(uri, diagnostics)

if __name__ == "__main__":
    server.start_io()
