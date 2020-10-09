
import AreaModification from "../recast/AreaModification.js"
class SampleAreaModifications {

static SAMPLE_POLYAREA_TYPE_MASK = 0x07;
static SAMPLE_POLYAREA_TYPE_GROUND = 0x1;
static SAMPLE_POLYAREA_TYPE_WATER = 0x2;
static SAMPLE_POLYAREA_TYPE_ROAD = 0x3;
static SAMPLE_POLYAREA_TYPE_DOOR = 0x4;
static SAMPLE_POLYAREA_TYPE_GRASS = 0x5;
static SAMPLE_POLYAREA_TYPE_JUMP = 0x6;

static SAMPLE_AREAMOD_GROUND = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GROUND,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
static SAMPLE_AREAMOD_WATER = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_WATER,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
static SAMPLE_AREAMOD_ROAD = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_ROAD,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
static SAMPLE_AREAMOD_GRASS = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_GRASS,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_MASK);
static SAMPLE_AREAMOD_DOOR = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_DOOR);
static SAMPLE_AREAMOD_JUMP = new AreaModification(SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP,
  SampleAreaModifications.SAMPLE_POLYAREA_TYPE_JUMP);

static  SAMPLE_POLYFLAGS_WALK = 0x01;	// Ability to walk (ground, grass, road)
static  SAMPLE_POLYFLAGS_SWIM = 0x02;   // Ability to swim (water).
static  SAMPLE_POLYFLAGS_DOOR = 0x04;   // Ability to move through doors.
static  SAMPLE_POLYFLAGS_JUMP = 0x08;   // Ability to jump.
static  SAMPLE_POLYFLAGS_DISABLED = 0x10; // Disabled polygon
static  SAMPLE_POLYFLAGS_ALL = 0xfff; // All abilities.
}

export default SampleAreaModifications;
