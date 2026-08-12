# Configuration

Copy `config.example.R` to `config.R` and edit the local paths before
running the analysis.

```r
file.copy(
  "config/config.example.R",
  "config/config.R"
)
```

`config/config.R` is excluded from version control because it may contain  
machine-specific paths. Do not add restricted-data locations, credentials,  
or personal information to `config.example.R`.
