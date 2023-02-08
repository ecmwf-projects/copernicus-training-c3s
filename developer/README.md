# Developer instructions

# Create a JUPYTER-BUILD environment with dependencies installed:
```
. ./developer/create-environment-for-build
```

## Rebuild jupyter book

```
rm -rf _build
jupyter-book build --all .
```

## Copy to gh-pages

```
ghp-import -n -p -f _build/html
```