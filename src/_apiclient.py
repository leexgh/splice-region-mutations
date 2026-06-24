"""Shared API client utilities.

Two features matter for re-runnability:
1. **Per-API JSON cache** at `out/cache/{api}/{key}.json`. Successful responses
   are written here verbatim; the next run reads from cache instead of calling
   the API. Cache key is caller-supplied (typically the 5-column variant_key).
2. **Polite throttle + exponential backoff on 429 / 5xx**. Default delay is
   conservative per API; override via `delay` argument.

Helpers:
- `variant_key(chrom, start, end, ref, alt) -> "1_12345_12345_A_T"` — canonical
  filename-safe per-variant cache key shared across all annotation steps.
- `polite_get(url, *, api, key, ...) -> dict | None` — cached GET that returns
  parsed JSON (or `None` for caller-tagged "not found" responses, e.g. 404).
- `polite_post(url, *, api, key, json_body, ...)` — same for POST.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Optional

import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CACHE_ROOT = PROJECT_ROOT / "out" / "cache"

log = logging.getLogger("apiclient")
if not log.handlers:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("[%(name)s] %(message)s"))
    log.addHandler(h)
    log.setLevel(logging.INFO)


import hashlib


def variant_key(chrom: str, start: int, end: int, ref: str, alt: str) -> str:
    """Filename-safe canonical key shared across the pipeline.

    For SNVs / short indels we keep the human-readable form. When ref or alt
    is long (e.g., a large deletion's reference sequence), the human-readable
    form would exceed the OS filename limit (255 bytes on macOS HFS+/APFS),
    so we substitute a 16-char SHA1 prefix for any allele >24 chars and
    record the full alleles inside the cached payload instead.
    """
    def shorten(s: str) -> str:
        if len(s) <= 24:
            return s
        return f"len{len(s)}_{hashlib.sha1(s.encode()).hexdigest()[:16]}"

    return f"{chrom}_{start}_{end}_{shorten(ref)}_{shorten(alt)}"


def _cache_path(api: str, key: str) -> Path:
    p = CACHE_ROOT / api / f"{key}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def cached_read(api: str, key: str) -> Optional[Any]:
    """Return cached payload or None."""
    p = _cache_path(api, key)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            log.warning("cache %s corrupt, will refetch", p)
            return None
    return None


def cached_write(api: str, key: str, payload: Any) -> None:
    _cache_path(api, key).write_text(json.dumps(payload, indent=None, separators=(",", ":")))


def _backoff_request(
    method: str,
    url: str,
    *,
    headers: Optional[dict] = None,
    params: Optional[dict] = None,
    json_body: Optional[Any] = None,
    timeout: float = 30.0,
    max_retries: int = 5,
    initial_backoff: float = 2.0,
) -> requests.Response:
    """One request with exponential backoff on 429 / 5xx."""
    sleep = initial_backoff
    last: Optional[requests.Response] = None
    for attempt in range(max_retries):
        last = requests.request(
            method,
            url,
            headers=headers,
            params=params,
            json=json_body,
            timeout=timeout,
        )
        if last.status_code < 500 and last.status_code != 429:
            return last
        log.warning(
            "  %s %s -> HTTP %s, retry in %.1fs (attempt %d/%d)",
            method,
            url[:80],
            last.status_code,
            sleep,
            attempt + 1,
            max_retries,
        )
        time.sleep(sleep)
        sleep *= 2
    assert last is not None
    return last


def polite_get(
    url: str,
    *,
    api: str,
    key: str,
    params: Optional[dict] = None,
    headers: Optional[dict] = None,
    delay: float = 0.5,
    not_found_codes: tuple[int, ...] = (404,),
    timeout: float = 30.0,
) -> Optional[Any]:
    """Cached, throttled GET. Returns parsed JSON, or None if response is a
    documented "not found" (e.g. 404). Raises on other 4xx after retries."""
    cached = cached_read(api, key)
    if cached is not None:
        return cached if cached != {"__not_found__": True} else None

    time.sleep(delay)
    resp = _backoff_request("GET", url, headers=headers, params=params, timeout=timeout)
    if resp.status_code in not_found_codes:
        cached_write(api, key, {"__not_found__": True})
        return None
    if not resp.ok:
        raise RuntimeError(f"{api} GET {url} -> {resp.status_code}: {resp.text[:200]}")
    payload = resp.json()
    cached_write(api, key, payload)
    return payload


def polite_post(
    url: str,
    *,
    api: str,
    key: str,
    json_body: Any,
    headers: Optional[dict] = None,
    delay: float = 0.5,
    not_found_codes: tuple[int, ...] = (404,),
    timeout: float = 60.0,
) -> Optional[Any]:
    """Cached, throttled POST."""
    cached = cached_read(api, key)
    if cached is not None:
        return cached if cached != {"__not_found__": True} else None

    time.sleep(delay)
    resp = _backoff_request(
        "POST",
        url,
        headers=headers,
        json_body=json_body,
        timeout=timeout,
    )
    if resp.status_code in not_found_codes:
        cached_write(api, key, {"__not_found__": True})
        return None
    if not resp.ok:
        raise RuntimeError(f"{api} POST {url} -> {resp.status_code}: {resp.text[:300]}")
    payload = resp.json()
    cached_write(api, key, payload)
    return payload
