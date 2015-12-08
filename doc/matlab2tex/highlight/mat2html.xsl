<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://www.w3.org/TR/xhtml1/strict">
<xsl:preserve-space elements="line" />
<xsl:output
   method="html"
   indent="yes" 
   encoding="iso-8859-1"/>
   
<xsl:template match="mfile">
 <html>
  <head>
   <title>
    <xsl:value-of select="@name"/>
   </title>
   <meta name="generator" content="highlight.m (c) 2003 Guillaume Flandin"/>
   <style type="text/css">
    .comment {color: #228B22;}
    .string {color: #B20000;}
    .keyword, .cont {color: #0000FF;}
    .cont {text-decoration: underline;}
    .code {color: #000000;}
   </style>
  </head>
  <body>
   <pre class="mcode">
    <xsl:apply-templates/>
   </pre>
  </body>
 </html>
</xsl:template>

<xsl:template match="line">
	<xsl:value-of select="@nb"/>
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="keyword">
  <span class="keyword">
    <xsl:apply-templates/>
  </span>
</xsl:template>

<xsl:template match="cont">
  <span class="cont">
    <xsl:apply-templates/>
  </span>
</xsl:template>

<xsl:template match="comment">
  <span class="comment">
    <xsl:apply-templates/>
  </span>
</xsl:template>

<xsl:template match="string">
  <span class="string">
    <xsl:apply-templates/>
  </span>
</xsl:template>

</xsl:stylesheet>

<!--
%**************************************************************************************
%*                         Stylesheet Author: Guillaume Flandin                       *
%*                                      June 2003                                     *
%*      PLEASE KEEP THIS NOTE WHEN USING, MODIFING OR AMELIORATING THIS SHEET         *
%*                   Send your comments to Guillaume@artefact.tk                      *
%**************************************************************************************
-->
