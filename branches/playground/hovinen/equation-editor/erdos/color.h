#ifndef BONOBO_AUDIO_ULAW_COLOR_H
#define BONOBO_AUDIO_ULAW_COLOR_H

void     color_init      (void);

/* Return the pixel value for the given red, green and blue */
int      color_alloc      (gushort red, gushort green, gushort blue);
void     color_alloc_name (const char *name, GdkColor *color);
void     color_alloc_gdk  (GdkColor *color);

/* Colors used by any audio-ulaw item */
extern GdkColor ia_white, ia_black;

#endif /* BONOBO_AUDIO_ULAW_COLOR_H */
