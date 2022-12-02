pandoc $1 \
    --from markdown \
    --standalone \
    --toc \
    --toc-depth 2 \
    --shift-heading-level-by=-1 \
    --to html \
    --template=./html_templates/uikit.html \
    --katex \
    --output $2 

