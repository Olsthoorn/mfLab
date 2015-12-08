<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output omit-xml-declaration="yes" />

<xsl:template match="mfile">
 \documentclass[a4paper,10pt]{article}
  \usepackage{alltt}
  \usepackage{color}
  \usepackage{fullpage}
  \definecolor{string}{rgb}{0.7,0.0,0.0}
  \definecolor{comment}{rgb}{0.13,0.54,0.13}
  \definecolor{keyword}{rgb}{0.0,0.0,1.0}
  \title{<xsl:value-of select="@name"/>}
  \author{\textsc{Matlab}, The Mathworks, Inc.}
  \begin{document}
  \maketitle
  \begin{alltt}
    <xsl:apply-templates/>
  \end{alltt}
 \end{document}
</xsl:template>

<xsl:template match="line">
  <xsl:value-of select="@nb"/><xsl:apply-templates/>
</xsl:template>

<xsl:template match="keyword">
  \textcolor{keyword}{<xsl:apply-templates/>}
</xsl:template>

<xsl:template match="comment">
  \textcolor{comment}{<xsl:apply-templates/>}
</xsl:template>

<xsl:template match="string">
  \textcolor{string}{<xsl:apply-templates/>}
</xsl:template>

</xsl:stylesheet>

<!--
There's missing a recursive function that generates LaTeX special characters.
It may be useful to glance at:
http://opera.inrialpes.fr/people/Tayeb.Lemlouma/MULTIMEDIA/XSLT/XML2LaTeX.xsl
-->
