rule ffmpeg_extra:
    """
    Created:
        2018-10-23 10:24:09
    Aim:
        Trim video without rendering (very fast)
    Note:
        https://vollnixx.wordpress.com/2012/06/01/howto-cut-a-video-directly-with-ffmpeg-without-transcoding/
        You only need the -ss and -t options. -ss is the starttime and -t the duration, so if you want 10 seconds from a video starting from minute one, you would use this.
        You can use seconds or hh:mm:ss[.xxx] as arguments for the start time and duration, i prefer the second option.
        The options -vcodec copy and -acodec copy are used to only cut the video and disable the re-encoding, this really speed things up.
    Test:
        out/ffmpeg/-vcodec_copy_-acodec_copy_-copyinkf_-ss_00:00:00.000_-t_00:00:30.000/inp/danse/salsa/tropicuba/1.mp4
    """
    input:
        video="out/{filler}.{ext}"
    output:
        video="out/ffmpeg/{extra}/{filler}.{ext}"
    params:
        extra=params_extra
    wildcard_constraints:
        ext="mp4|webm" #maybe more
    conda:
        "../envs/ffmpeg.yaml"
    shell:
        """
        ffmpeg -i {input.video} {params.extra} {output.video}
        """

rule ffmpeg_vcodec_acodec_ss_t:
    """
    Created:
        2017-09-27 09:48:40
    Aim:
        Trim video without rendering (very fast)
    Note:
        https://vollnixx.wordpress.com/2012/06/01/howto-cut-a-video-directly-with-ffmpeg-without-transcoding/
        You only need the -ss and -t options. -ss is the starttime and -t the duration, so if you want 10 seconds from a video starting from minute one, you would use this.
        You can use seconds or hh:mm:ss[.xxx] as arguments for the start time and duration, i prefer the second option.
        The options -vcodec copy and -acodec copy are used to only cut the video and disable the re-encoding, this really speed things up.
    Test:
        out/ffmpeg/vcodec-copy_acodec-copy_ss-00:00:00.000_t-00:03:47.900/inp/dwhelper/Evo-Devo-Despacito_Biology_Parody-A_Capella_Science.webm
        out/ffmpeg/vcodec-copy_acodec-copy_ss-00:00:00.000_t-00:00:30.000/inp/danse/salsa/tropicuba/1.mp4

    """
    input:
        video="out/{filler}.{ext}"
    output:
        video="out/ffmpeg/vcodec-{vcodec}_acodec-{acodec}_ss-{ss}_t-{t}/{filler}.{ext}"
    wildcard_constraints:
        vcodec="copy", #look doc someday to add other arguments.
        acodec="copy",
        ss="[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}",
        t="[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}",
        ext="mp4|webm" #maybe more
    conda:
        "../envs/ffmpeg.yaml"
    shell:
        """
        ffmpeg -i {input.video} \
            -vcodec {wildcards.vcodec} \
            -acodec {wildcards.acodec} \
            -ss {wildcards.ss} \
            -t {wildcards.t} \
            {output.video}
        """
    
