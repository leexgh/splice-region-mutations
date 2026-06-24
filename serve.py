#!/usr/bin/env python3
"""Start a local HTTP server to serve the splice-region report with IGV.js support.

IGV.js requires HTTP (not file://) to load reference genomes from CDN.
Usage: python serve.py [port]
"""
import http.server
import os
import sys
import threading
import webbrowser

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

os.chdir(ROOT)

handler = http.server.SimpleHTTPRequestHandler

def open_browser():
    webbrowser.open(f"http://localhost:{PORT}/report.html")

threading.Timer(0.5, open_browser).start()
print(f"Serving http://localhost:{PORT}/report.html")
print("Press Ctrl+C to stop.")
http.server.HTTPServer(("", PORT), handler).serve_forever()
