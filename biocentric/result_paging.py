"""Bound the number of result rows sent to each browser table."""


def page_rows(st, rows, key, label, page_size=100):
    if len(rows) <= page_size:
        return rows
    pages = (len(rows) + page_size - 1) // page_size
    page = st.selectbox(f"{label} page", list(range(pages)),
                        format_func=lambda value: f"{value + 1} of {pages}", key=key)
    start = page * page_size
    st.caption(f"Showing {start + 1}–{min(start + page_size, len(rows))} of {len(rows)} records. Downloads include every recorded row.")
    return rows[start:start + page_size]
