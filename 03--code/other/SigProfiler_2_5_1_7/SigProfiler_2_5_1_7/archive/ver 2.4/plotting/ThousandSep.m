function out = ThousandSep(in)
    import java.text.*
    v = DecimalFormat;
    out = char(v.format(in));
end